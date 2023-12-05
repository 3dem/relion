#ifndef BP_GPU_KERNELS_H_
#define BP_GPU_KERNELS_H_

#include <cassert>
#include <sycl/sycl.hpp>

#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/sycl/sycl_settings.h"
#include "src/acc/sycl/sycl_kernels/helper.h"
#include "src/acc/sycl/sycl_kernels/helper_gpu.h"

namespace syclGpuKernels
{

template <bool CTF_PREMULTIPLIED, bool SGD>
void sycl_kernel_backproject2D(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		int trans_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers, XFLOAT *g_model_real,
		XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor,
		int img_x, int img_y, int img_xy, int mdl_x, int mdl_inity,
		XFLOAT *s_eulers)
{
	const int tid = nit.get_local_id(2);
	const int img = nit.get_group_linear_id();

	const int img_y_half = img_y / 2;
	const XFLOAT max_r2_out = sycl::floor(max_r2 * padding_factor * padding_factor);
	const XFLOAT inv_weight_norm = 1.0f / weight_norm;

	if (tid == 0)
		s_eulers[0] = g_eulers[img*9 + 0];
	else if (tid == 1)
		s_eulers[1] = g_eulers[img*9 + 1];
	else if (tid == 2)
		s_eulers[2] = g_eulers[img*9 + 3];
	else if (tid == 3)
		s_eulers[3] = g_eulers[img*9 + 4];

	__group_barrier(nit);

	const int pixel_pass_num = img_xy/BP_2D_BLOCK_SIZE + 1;
	for (int pass = 0; pass < pixel_pass_num; pass++)
	{
		const int pixel = pass*BP_2D_BLOCK_SIZE + tid;

		if (pixel >= img_xy)
			continue;

		int x = pixel % img_x;
		int y = pixel / img_x;
		if (y > img_y_half)
			y -= img_y;

		XFLOAT minvsigma2 = g_Minvsigma2s[pixel];
		XFLOAT img_real = g_img_real[pixel];
		XFLOAT img_imag = g_img_imag[pixel];
		XFLOAT ctf = g_ctfs[pixel];

		XFLOAT ref_real = 0.0f;
		XFLOAT ref_imag = 0.0f;
		if (SGD)
		{
			projector.project2Dmodel(x, y, s_eulers[0], s_eulers[1], s_eulers[2], s_eulers[3], ref_real, ref_imag);
			ref_real *= ctf;
			ref_imag *= ctf;
		}

		XFLOAT Fweight = 0.0f;
		XFLOAT real = 0.0f;
		XFLOAT imag = 0.0f;
		for (int itrans = 0; itrans < trans_num; itrans++)
		{
			XFLOAT weight = g_weights[img*trans_num + itrans];

			if (weight >= significant_weight)
			{
				if (CTF_PREMULTIPLIED)
					weight = weight * inv_weight_norm * minvsigma2;
				else
					weight = weight * inv_weight_norm * ctf * minvsigma2;

				Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!

				XFLOAT temp_real, temp_imag;
				translatePixel(x, y, g_trans_x[itrans], g_trans_y[itrans], img_real, img_imag, temp_real, temp_imag);

				if (SGD)
				{
					real += (temp_real - ref_real) * weight;
					imag += (temp_imag - ref_imag) * weight;
				}
				else
				{
					real += temp_real * weight;
					imag += temp_imag * weight;
				}
			}
		}

		if (Fweight > 0.0f)
		{
			// Get logical coordinates in the 3D map
			XFLOAT xp = (s_eulers[0]*x + s_eulers[1]*y) * padding_factor;
			XFLOAT yp = (s_eulers[2]*x + s_eulers[3]*y) * padding_factor;

			// Only consider pixels that are projected inside the allowed circle in output coordinates.
			//     --JZ, Nov. 26th 2018
			if ((xp*xp + yp*yp) > max_r2_out)	// TODO: This is very subtle and make big difference for some cases
				continue;

			// Only asymmetric half is stored
			if (xp < 0.0f)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				imag = -imag;
			}

			int x0 = sycl::floor(xp);
			XFLOAT fx = xp - x0;

			int y0 = sycl::floor(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;

			int offset1 = y0*mdl_x + x0;
			int offset2 = offset1 + 1;
			int offset3 = offset1 + mdl_x;
			int offset4 = offset1 + mdl_x + 1;

			XFLOAT mfx = 1.0f - fx;
			XFLOAT mfy = 1.0f - fy;

			XFLOAT dd00 = mfy * mfx;
			XFLOAT dd01 = mfy *  fx;
			XFLOAT dd10 =  fy * mfx;
			XFLOAT dd11 =  fy *  fx;

			atomic_ref_dev( g_model_real  [offset1] ).fetch_add(dd00 * real);
			atomic_ref_dev( g_model_imag  [offset1] ).fetch_add(dd00 * imag);
			atomic_ref_dev( g_model_weight[offset1] ).fetch_add(dd00 * Fweight);

			atomic_ref_dev( g_model_real  [offset2] ).fetch_add(dd01 * real);
			atomic_ref_dev( g_model_imag  [offset2] ).fetch_add(dd01 * imag);
			atomic_ref_dev( g_model_weight[offset2] ).fetch_add(dd01 * Fweight);

			atomic_ref_dev( g_model_real  [offset3] ).fetch_add(dd10 * real);
			atomic_ref_dev( g_model_imag  [offset3] ).fetch_add(dd10 * imag);
			atomic_ref_dev( g_model_weight[offset3] ).fetch_add(dd10 * Fweight);

			atomic_ref_dev( g_model_real  [offset4] ).fetch_add(dd11 * real);
			atomic_ref_dev( g_model_imag  [offset4] ).fetch_add(dd11 * imag);
			atomic_ref_dev( g_model_weight[offset4] ).fetch_add(dd11 * Fweight);
		}
	}
}

template <bool DATA3D, bool CTF_PREMULTIPLIED, bool SGD>
void sycl_kernel_backproject3D(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		int trans_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real, XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor, int img_x,
		int img_y, int img_z, int img_xyz, int mdl_x,
		int mdl_y, int mdl_inity, int mdl_initz,
		XFLOAT *s_eulers)
{
	const int tid = nit.get_local_id(2);
	const int img = nit.get_group_linear_id();

	const int img_y_half = img_y / 2;
	const int img_z_half = img_z / 2;
	const XFLOAT max_r2_vol = sycl::floor(max_r2 * padding_factor * padding_factor);
	const XFLOAT inv_weight_norm = 1.0f / weight_norm;

	if (tid < 9)
		s_eulers[tid] = g_eulers[img*9 + tid];

	__group_barrier(nit);

	int pixel_pass_num;
	if (DATA3D)
		pixel_pass_num = img_xyz/BP_DATA3D_BLOCK_SIZE + 1;
	else
		pixel_pass_num = img_xyz/BP_REF3D_BLOCK_SIZE + 1;

	for (int pass = 0; pass < pixel_pass_num; pass++)
	{
		int pixel;
		if (DATA3D)
			pixel = pass*BP_DATA3D_BLOCK_SIZE + tid;
		else
			pixel = pass*BP_REF3D_BLOCK_SIZE + tid;

		if (pixel >= img_xyz)
			continue;

		int x, y, z, xy;
		if (DATA3D)
		{
			z =  pixel / (img_x*img_y);
			xy = pixel % (img_x*img_y);
			x =     xy % img_x;
			y =     xy / img_x;

			if (z > img_z_half)
			{
				z = z - img_z;

				if (x == 0)
					continue;
			}
		}
		else
		{
			x = pixel % img_x;
			y = pixel / img_x;
		}
		if (y > img_y_half)
			y = y - img_y;

		//WAVG
		XFLOAT minvsigma2 = g_Minvsigma2s[pixel];
		XFLOAT img_real = g_img_real[pixel];
		XFLOAT img_imag = g_img_imag[pixel];
		XFLOAT ctf = g_ctfs[pixel];

		XFLOAT ref_real = 0.0f;
		XFLOAT ref_imag = 0.0f;
		if (SGD)
		{
			if (DATA3D)
				projector.project3Dmodel(
					x, y, z,
					s_eulers[0], s_eulers[1], s_eulers[2],
					s_eulers[3], s_eulers[4], s_eulers[5],
					s_eulers[6], s_eulers[7], s_eulers[8],
					ref_real, ref_imag);
			else
				projector.project3Dmodel(
					x, y,
					s_eulers[0], s_eulers[1],
					s_eulers[3], s_eulers[4],
					s_eulers[6], s_eulers[7],
					ref_real, ref_imag);

			ref_real *= ctf;
			ref_imag *= ctf;
		}

		XFLOAT Fweight = 0.0f;
		XFLOAT real = 0.0f;
		XFLOAT imag = 0.0f;
		for (int itrans = 0; itrans < trans_num; itrans++)
		{
			XFLOAT weight = g_weights[img*trans_num + itrans];

			if (weight >= significant_weight)
			{
				if (CTF_PREMULTIPLIED)
					weight = weight * inv_weight_norm * minvsigma2;
				else
					weight = weight * inv_weight_norm * ctf * minvsigma2;

				Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!

				XFLOAT temp_real, temp_imag;
				if (DATA3D)
					translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
				else
					translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

				if (SGD)
				{
					real += (temp_real - ref_real) * weight;
					imag += (temp_imag - ref_imag) * weight;
				}
				else
				{
					real += temp_real * weight;
					imag += temp_imag * weight;
				}
			}
		}

		//BP
		if (Fweight > 0.0f)
		{
			// Get logical coordinates in the 3D map
			XFLOAT xp, yp, zp;
			if (DATA3D)
			{
				xp = (s_eulers[0]*x + s_eulers[1]*y + s_eulers[2]*z) * padding_factor;
				yp = (s_eulers[3]*x + s_eulers[4]*y + s_eulers[5]*z) * padding_factor;
				zp = (s_eulers[6]*x + s_eulers[7]*y + s_eulers[8]*z) * padding_factor;
			}
			else
			{
				xp = (s_eulers[0]*x + s_eulers[1]*y) * padding_factor;
				yp = (s_eulers[3]*x + s_eulers[4]*y) * padding_factor;
				zp = (s_eulers[6]*x + s_eulers[7]*y) * padding_factor;
			}

			// Only consider pixels that are projected inside the sphere in output coordinates.
			//     --JZ, Oct. 18. 2018
			if ((xp*xp + yp*yp + zp*zp) > max_r2_vol)	// TODO: This is very subtle and make big difference for some cases
				continue;

			// Only asymmetric half is stored
			if (xp < 0.0f)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				imag = -imag;
			}

			int x0 = sycl::floor(xp);
			XFLOAT fx = xp - x0;

			int y0 = sycl::floor(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;

			int z0 = sycl::floor(zp);
			XFLOAT fz = zp - z0;
			z0 -= mdl_initz;

			int offset1 = z0*mdl_x*mdl_y + y0*mdl_x + x0;
			int offset2 = offset1 + 1;
			int offset3 = offset1 + mdl_x;
			int offset4 = offset1 + mdl_x + 1;
			int offset5 = offset1 + mdl_x*mdl_y;
			int offset6 = offset1 + mdl_x*mdl_y + 1;
			int offset7 = offset1 + mdl_x*mdl_y + mdl_x;
			int offset8 = offset1 + mdl_x*mdl_y + mdl_x + 1;

			XFLOAT mfx = 1.0f - fx;
			XFLOAT mfy = 1.0f - fy;
			XFLOAT mfz = 1.0f - fz;

			XFLOAT dd000 = mfz * mfy * mfx;

			atomic_ref_dev( g_model_real  [offset1] ).fetch_add(dd000 * real);
			atomic_ref_dev( g_model_imag  [offset1] ).fetch_add(dd000 * imag);
			atomic_ref_dev( g_model_weight[offset1] ).fetch_add(dd000 * Fweight);

			XFLOAT dd001 = mfz * mfy *  fx;

			atomic_ref_dev( g_model_real  [offset2] ).fetch_add(dd001 * real);
			atomic_ref_dev( g_model_imag  [offset2] ).fetch_add(dd001 * imag);
			atomic_ref_dev( g_model_weight[offset2] ).fetch_add(dd001 * Fweight);

			XFLOAT dd010 = mfz *  fy * mfx;

			atomic_ref_dev( g_model_real  [offset3] ).fetch_add(dd010 * real);
			atomic_ref_dev( g_model_imag  [offset3] ).fetch_add(dd010 * imag);
			atomic_ref_dev( g_model_weight[offset3] ).fetch_add(dd010 * Fweight);

			XFLOAT dd011 = mfz *  fy *  fx;

			atomic_ref_dev( g_model_real  [offset4] ).fetch_add(dd011 * real);
			atomic_ref_dev( g_model_imag  [offset4] ).fetch_add(dd011 * imag);
			atomic_ref_dev( g_model_weight[offset4] ).fetch_add(dd011 * Fweight);

			XFLOAT dd100 =  fz * mfy * mfx;

			atomic_ref_dev( g_model_real  [offset5] ).fetch_add(dd100 * real);
			atomic_ref_dev( g_model_imag  [offset5] ).fetch_add(dd100 * imag);
			atomic_ref_dev( g_model_weight[offset5] ).fetch_add(dd100 * Fweight);

			XFLOAT dd101 =  fz * mfy *  fx;

			atomic_ref_dev( g_model_real  [offset6] ).fetch_add(dd101 * real);
			atomic_ref_dev( g_model_imag  [offset6] ).fetch_add(dd101 * imag);
			atomic_ref_dev( g_model_weight[offset6] ).fetch_add(dd101 * Fweight);

			XFLOAT dd110 =  fz *  fy * mfx;

			atomic_ref_dev( g_model_real  [offset7] ).fetch_add(dd110 * real);
			atomic_ref_dev( g_model_imag  [offset7] ).fetch_add(dd110 * imag);
			atomic_ref_dev( g_model_weight[offset7] ).fetch_add(dd110 * Fweight);

			XFLOAT dd111 =  fz *  fy *  fx;

			atomic_ref_dev( g_model_real  [offset8] ).fetch_add(dd111 * real);
			atomic_ref_dev( g_model_imag  [offset8] ).fetch_add(dd111 * imag);
			atomic_ref_dev( g_model_weight[offset8] ).fetch_add(dd111 * Fweight);
		}
	}
}

} // namespace

#endif
