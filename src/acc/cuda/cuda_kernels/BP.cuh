#ifndef CUDA_BP_KERNELS_CUH_
#define CUDA_BP_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"


/*
 *   	BP KERNELS
 */

template < bool CTF_PREMULTIPLIED >
__global__ void cuda_kernel_backproject2D(
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT* g_weights,
		XFLOAT* g_Minvsigma2s,
		XFLOAT* g_ctfs,
		unsigned long translation_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_xy,
		unsigned mdl_x,
		int mdl_inity)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;

	int img_y_half = img_y / 2;
	int max_r2_out = max_r2 * padding_factor * padding_factor;

	__shared__ XFLOAT s_eulers[4];

	XFLOAT minvsigma2, ctf, img_real, img_imag, Fweight, real, imag, weight;

	if (tid == 0)
		s_eulers[0] = g_eulers[img*9+0] * padding_factor;
	else if (tid == 1)
		s_eulers[1] = g_eulers[img*9+1] * padding_factor;
	else if (tid == 2)
		s_eulers[2] = g_eulers[img*9+3] * padding_factor;
	else if (tid == 3)
		s_eulers[3] = g_eulers[img*9+4] * padding_factor;

	__syncthreads();

	int pixel_pass_num(ceilf((float)img_xy/(float)BP_2D_BLOCK_SIZE));

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{
		unsigned pixel = (pass * BP_2D_BLOCK_SIZE) + tid;

		if (pixel >= img_xy)
			continue;

		int x = pixel % img_x;
		int y = (int)floorf( (float)pixel / (float)img_x);

		if (y > img_y_half)
		{
			y -= img_y;
		}

		//WAVG
		minvsigma2 = __ldg(&g_Minvsigma2s[pixel]);
		ctf = __ldg(&g_ctfs[pixel]);
		img_real = __ldg(&g_img_real[pixel]);
		img_imag = __ldg(&g_img_imag[pixel]);
		Fweight = (XFLOAT) 0.0;
		real = (XFLOAT) 0.0;
		imag = (XFLOAT) 0.0;

		XFLOAT temp_real, temp_imag;

		for (unsigned long itrans = 0; itrans < translation_num; itrans++)
		{
			weight = g_weights[img * translation_num + itrans];

			if (weight >= significant_weight)
			{
				if(CTF_PREMULTIPLIED)
				{
					weight = (weight / weight_norm) * minvsigma2;
					Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!
				}
				else
				{
					weight = (weight / weight_norm) * ctf * minvsigma2;
					Fweight += weight * ctf;
				}

				translatePixel(x, y, g_trans_x[itrans], g_trans_y[itrans], img_real, img_imag, temp_real, temp_imag);

				real += temp_real * weight;
				imag += temp_imag * weight;

			}
		}

		if (Fweight > (XFLOAT) 0.0)
		{

			// Get logical coordinates in the 3D map
			XFLOAT xp = (s_eulers[0] * x + s_eulers[1] * y );
			XFLOAT yp = (s_eulers[2] * x + s_eulers[3] * y );

			// Only consider pixels that are projected inside the allowed circle in output coordinates.
			//     --JZ, Nov. 26th 2018
			if ( ( xp * xp + yp * yp ) > max_r2_out )
				continue;

			// Only asymmetric half is stored
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			XFLOAT mfx = (XFLOAT) 1.0 - fx;
			XFLOAT mfy = (XFLOAT) 1.0 - fy;

			XFLOAT dd00 = mfy * mfx;
			XFLOAT dd01 = mfy *  fx;
			XFLOAT dd10 =  fy * mfx;
			XFLOAT dd11 =  fy *  fx;

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x0], dd00 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x0], dd00 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x0], dd00 * Fweight);

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x1], dd01 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x1], dd01 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x1], dd01 * Fweight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x0], dd10 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x0], dd10 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x0], dd10 * Fweight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x1], dd11 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x1], dd11 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x1], dd11 * Fweight);
		}
	}
}

template < bool DATA3D, bool CTF_PREMULTIPLIED >
__global__ void cuda_kernel_backproject3D(
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		XFLOAT* g_weights,
		XFLOAT* g_Minvsigma2s,
		XFLOAT* g_ctfs,
		unsigned long translation_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_z,
		unsigned img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;
	
	int img_y_half = img_y / 2;
	int img_z_half = img_z / 2;
	
	int max_r2_vol = max_r2 * padding_factor * padding_factor;

	__shared__ XFLOAT s_eulers[9];
	XFLOAT minvsigma2, ctf, img_real, img_imag, Fweight, real, imag, weight;

	if (tid < 9)
		s_eulers[tid] = g_eulers[img*9+tid];

	__syncthreads();

	int pixel_pass_num(0);
	if(DATA3D)
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_DATA3D_BLOCK_SIZE));
	else
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_REF3D_BLOCK_SIZE));

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
    {
		unsigned pixel(0);
		if(DATA3D)
			pixel = (pass * BP_DATA3D_BLOCK_SIZE) + tid;
		else
			pixel = (pass * BP_REF3D_BLOCK_SIZE) + tid;

		if (pixel >= img_xyz)
			continue;

		int x,y,z,xy;

		if(DATA3D)
		{
			z =  floorfracf(pixel, img_x*img_y);
			xy = pixel % (img_x*img_y);
			x =             xy  % img_x;
			y = floorfracf( xy,   img_x);
			
			if (z > img_z_half)
			{
				z = z - img_z;

				if(x==0)
					continue;
			}
		}
		else
		{
			x =             pixel % img_x;
			y = floorfracf( pixel , img_x);
		}
		if (y > img_y_half)
		{
			y = y - img_y;
		}
		
		//WAVG
		minvsigma2 = __ldg(&g_Minvsigma2s[pixel]);
		ctf = __ldg(&g_ctfs[pixel]);
		img_real = __ldg(&g_img_real[pixel]);
		img_imag = __ldg(&g_img_imag[pixel]);
		Fweight = (XFLOAT) 0.0;
		real = (XFLOAT) 0.0;
		imag = (XFLOAT) 0.0;

		XFLOAT temp_real, temp_imag;

		for (unsigned long itrans = 0; itrans < translation_num; itrans++)
		{
			weight = g_weights[img * translation_num + itrans];

			if (weight >= significant_weight)
			{
				if(CTF_PREMULTIPLIED)
				{
					weight = (weight / weight_norm) * minvsigma2;
					Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!
				}
				else
				{
					weight = (weight / weight_norm) * ctf * minvsigma2;
					Fweight += weight * ctf;
				}
				
				if(DATA3D)
					translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
				else
					translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

				real += temp_real * weight;
				imag += temp_imag * weight;
			}
		}

		//BP
		if (Fweight > (XFLOAT) 0.0)
		{
			// Get logical coordinates in the 3D map

			XFLOAT xp,yp,zp;
			if(DATA3D)
			{
				xp = (s_eulers[0] * x + s_eulers[1] * y + s_eulers[2] * z) * padding_factor;
				yp = (s_eulers[3] * x + s_eulers[4] * y + s_eulers[5] * z) * padding_factor;
				zp = (s_eulers[6] * x + s_eulers[7] * y + s_eulers[8] * z) * padding_factor;
			}
			else
			{
				xp = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
				yp = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
				zp = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;
			}
			
			// Only consider pixels that are projected inside the sphere in output coordinates.
			//     --JZ, Oct. 18. 2018			
			if ( ( xp * xp + yp * yp  + zp * zp ) > max_r2_vol)
				continue;
		
			// Only asymmetric half is stored
			if (xp < (XFLOAT) 0.0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			int z0 = floorf(zp);
			XFLOAT fz = zp - z0;
			z0 -= mdl_initz;
			int z1 = z0 + 1;

			XFLOAT mfx = (XFLOAT)1.0 - fx;
			XFLOAT mfy = (XFLOAT)1.0 - fy;
			XFLOAT mfz = (XFLOAT)1.0 - fz;

			XFLOAT dd000 = mfz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * Fweight);

			XFLOAT dd001 = mfz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * Fweight);

			XFLOAT dd010 = mfz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * Fweight);

			XFLOAT dd011 = mfz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * Fweight);

			XFLOAT dd100 =  fz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * Fweight);

			XFLOAT dd101 =  fz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * Fweight);

			XFLOAT dd110 =  fz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * Fweight);

			XFLOAT dd111 =  fz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * Fweight);

		}
	}
}


template < bool DATA3D, bool CTF_PREMULTIPLIED >
__global__ void cuda_kernel_backproject3D_SGD(
		AccProjectorKernel projector,
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		XFLOAT* g_weights,
		XFLOAT* g_Minvsigma2s,
		XFLOAT* g_ctfs,
		unsigned long translation_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_z,
		unsigned img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;
	
	int img_y_half = img_y / 2;
	int img_z_half = img_z / 2;
	
	int max_r2_vol = max_r2 * padding_factor * padding_factor;

	__shared__ XFLOAT s_eulers[9];
	XFLOAT minvsigma2, ctf, img_real, img_imag, Fweight, real, imag, weight;

	if (tid < 9)
		s_eulers[tid] = g_eulers[img*9+tid];

	__syncthreads();

	int pixel_pass_num(0);
	if(DATA3D)
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_DATA3D_BLOCK_SIZE));
	else
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_REF3D_BLOCK_SIZE));

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
    {
		unsigned pixel(0);
		if(DATA3D)
			pixel = (pass * BP_DATA3D_BLOCK_SIZE) + tid;
		else
			pixel = (pass * BP_REF3D_BLOCK_SIZE) + tid;

		if (pixel >= img_xyz)
			continue;

		int x,y,z,xy;

		if(DATA3D)
		{
			z =  floorfracf(pixel, img_x*img_y);
			xy = pixel % (img_x*img_y);
			x =             xy  % img_x;
			y = floorfracf( xy,   img_x);
			
			if (z > img_z_half)
			{
				z = z - img_z;

				if(x==0)
					continue;
			}
		}
		else
		{
			x =             pixel % img_x;
			y = floorfracf( pixel , img_x);
		}
		if (y > img_y_half)
		{
			y = y - img_y;
		}

		XFLOAT ref_real = (XFLOAT) 0.0;
		XFLOAT ref_imag = (XFLOAT) 0.0;

		if(DATA3D)
			projector.project3Dmodel(
				x,y,z,
				s_eulers[0], s_eulers[1], s_eulers[2],
				s_eulers[3], s_eulers[4], s_eulers[5],
				s_eulers[6], s_eulers[7], s_eulers[8],
				ref_real, ref_imag);
		else
			projector.project3Dmodel(
				x,y,
				s_eulers[0], s_eulers[1],
				s_eulers[3], s_eulers[4],
				s_eulers[6], s_eulers[7],
				ref_real, ref_imag);

		//WAVG
		minvsigma2 = __ldg(&g_Minvsigma2s[pixel]);
		ctf = __ldg(&g_ctfs[pixel]);
		img_real = __ldg(&g_img_real[pixel]);
		img_imag = __ldg(&g_img_imag[pixel]);
		Fweight = (XFLOAT) 0.0;
		real = (XFLOAT) 0.0;
		imag = (XFLOAT) 0.0;
		ref_real *= ctf;
		ref_imag *= ctf;

		XFLOAT temp_real, temp_imag;

		for (unsigned long itrans = 0; itrans < translation_num; itrans++)
		{
			weight = g_weights[img * translation_num + itrans];

			if (weight >= significant_weight)
			{
				if(CTF_PREMULTIPLIED)
				{
					weight = (weight / weight_norm) * minvsigma2;
					Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!
				}
				else
				{
					weight = (weight / weight_norm) * ctf * minvsigma2;
					Fweight += weight * ctf;
				}
				if(DATA3D)
					translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
				else
					translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

				real += (temp_real-ref_real) * weight;
				imag += (temp_imag-ref_imag) * weight;
			}
		}

		//BP
		if (Fweight > (XFLOAT) 0.0)
		{
			// Get logical coordinates in the 3D map

			XFLOAT xp,yp,zp;
			if(DATA3D)
			{
				xp = (s_eulers[0] * x + s_eulers[1] * y + s_eulers[2] * z) * padding_factor;
				yp = (s_eulers[3] * x + s_eulers[4] * y + s_eulers[5] * z) * padding_factor;
				zp = (s_eulers[6] * x + s_eulers[7] * y + s_eulers[8] * z) * padding_factor;
			}
			else
			{
				xp = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
				yp = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
				zp = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;
			}
			
			// Only consider pixels that are projected inside the sphere in output coordinates.
			//     --JZ, Nov. 26th 2018			
			if ( ( xp * xp + yp * yp  + zp * zp ) > max_r2_vol)
				continue;
			
			// Only asymmetric half is stored
			if (xp < (XFLOAT) 0.0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				zp = -zp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			int z0 = floorf(zp);
			XFLOAT fz = zp - z0;
			z0 -= mdl_initz;
			int z1 = z0 + 1;

			XFLOAT mfx = (XFLOAT)1.0 - fx;
			XFLOAT mfy = (XFLOAT)1.0 - fy;
			XFLOAT mfz = (XFLOAT)1.0 - fz;

			XFLOAT dd000 = mfz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x0], dd000 * Fweight);

			XFLOAT dd001 = mfz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x1], dd001 * Fweight);

			XFLOAT dd010 = mfz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x0], dd010 * Fweight);

			XFLOAT dd011 = mfz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * real);
			cuda_atomic_add(&g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * imag);
			cuda_atomic_add(&g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x1], dd011 * Fweight);

			XFLOAT dd100 =  fz * mfy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x0], dd100 * Fweight);

			XFLOAT dd101 =  fz * mfy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x1], dd101 * Fweight);

			XFLOAT dd110 =  fz *  fy * mfx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x0], dd110 * Fweight);

			XFLOAT dd111 =  fz *  fy *  fx;

			cuda_atomic_add(&g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * real);
			cuda_atomic_add(&g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * imag);
			cuda_atomic_add(&g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x1], dd111 * Fweight);

		}
	}
}


template < bool CTF_PREMULTIPLIED >
__global__ void cuda_kernel_backproject2D_SGD(
		AccProjectorKernel projector,
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT* g_weights,
		XFLOAT* g_Minvsigma2s,
		XFLOAT* g_ctfs,
		unsigned long translation_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		int max_r,
		int max_r2,
		XFLOAT padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_xy,
		unsigned mdl_x,
		int mdl_inity)
{
	unsigned tid = threadIdx.x;
	unsigned img = blockIdx.x;

	int img_y_half = img_y / 2;
	int max_r2_out = max_r2 * padding_factor * padding_factor;

	__shared__ XFLOAT s_eulers[4];

	XFLOAT minvsigma2, ctf, img_real, img_imag, Fweight, real, imag, weight;

	if (tid == 0)
		s_eulers[0] = g_eulers[img*9+0];
	else if (tid == 1)
		s_eulers[1] = g_eulers[img*9+1];
	else if (tid == 2)
		s_eulers[2] = g_eulers[img*9+3];
	else if (tid == 3)
		s_eulers[3] = g_eulers[img*9+4];

	__syncthreads();

	int pixel_pass_num(ceilf((float)img_xy/(float)BP_2D_BLOCK_SIZE));

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{
		unsigned pixel = (pass * BP_2D_BLOCK_SIZE) + tid;

		if (pixel >= img_xy)
			continue;

		int x = pixel % img_x;
		int y = (int)floorf( (float)pixel / (float)img_x);

		if (y > img_y_half)
		{
			y -= img_y;
		}

		//WAVG
		minvsigma2 = __ldg(&g_Minvsigma2s[pixel]);
		ctf = __ldg(&g_ctfs[pixel]);
		img_real = __ldg(&g_img_real[pixel]);
		img_imag = __ldg(&g_img_imag[pixel]);
		Fweight = (XFLOAT) 0.0;
		real = (XFLOAT) 0.0;
		imag = (XFLOAT) 0.0;

		XFLOAT temp_real, temp_imag;

		XFLOAT ref_real = (XFLOAT) 0.0;
		XFLOAT ref_imag = (XFLOAT) 0.0;

		projector.project2Dmodel(
			x,y,
			s_eulers[0], s_eulers[1],
			s_eulers[2], s_eulers[3],
			ref_real, ref_imag
		);
		ref_real *= ctf;
		ref_imag *= ctf;

		for (unsigned long itrans = 0; itrans < translation_num; itrans++)
		{
			weight = g_weights[img * translation_num + itrans];

			if (weight >= significant_weight)
			{
				if(CTF_PREMULTIPLIED)
				{
					weight = (weight / weight_norm) * minvsigma2;
					Fweight += weight * ctf; // SHWS 13feb2020: from now on when ctf_premultiplied, the ctf array actually contains ctf^2!
				}
				else
				{
					weight = (weight / weight_norm) * ctf * minvsigma2;
					Fweight += weight * ctf;
				}

				translatePixel(x, y, g_trans_x[itrans], g_trans_y[itrans], img_real, img_imag, temp_real, temp_imag);

				real += (temp_real-ref_real) * weight;
				imag += (temp_imag-ref_imag) * weight;
			}
		}

		if (Fweight > (XFLOAT) 0.0)
		{

			// Get logical coordinates in the 3D map
			XFLOAT xp = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
			XFLOAT yp = (s_eulers[2] * x + s_eulers[3] * y ) * padding_factor;

			// Only consider pixels that are projected inside the allowed circle in output coordinates.
			//     --JZ, Nov. 26th 2018
			if ( ( xp * xp + yp * yp ) > max_r2_out )
				continue;

			// Only asymmetric half is stored
			if (xp < 0)
			{
				// Get complex conjugated hermitian symmetry pair
				xp = -xp;
				yp = -yp;
				imag = -imag;
			}

			int x0 = floorf(xp);
			XFLOAT fx = xp - x0;
			int x1 = x0 + 1;

			int y0 = floorf(yp);
			XFLOAT fy = yp - y0;
			y0 -= mdl_inity;
			int y1 = y0 + 1;

			XFLOAT mfx = (XFLOAT) 1.0 - fx;
			XFLOAT mfy = (XFLOAT) 1.0 - fy;

			XFLOAT dd00 = mfy * mfx;
			XFLOAT dd01 = mfy *  fx;
			XFLOAT dd10 =  fy * mfx;
			XFLOAT dd11 =  fy *  fx;

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x0], dd00 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x0], dd00 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x0], dd00 * Fweight);

			cuda_atomic_add(&g_model_real  [y0 * mdl_x + x1], dd01 * real);
			cuda_atomic_add(&g_model_imag  [y0 * mdl_x + x1], dd01 * imag);
			cuda_atomic_add(&g_model_weight[y0 * mdl_x + x1], dd01 * Fweight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x0], dd10 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x0], dd10 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x0], dd10 * Fweight);

			cuda_atomic_add(&g_model_real  [y1 * mdl_x + x1], dd11 * real);
			cuda_atomic_add(&g_model_imag  [y1 * mdl_x + x1], dd11 * imag);
			cuda_atomic_add(&g_model_weight[y1 * mdl_x + x1], dd11 * Fweight);
		}
	}
}

#endif /* CUDA_PB_KERNELS_CUH_ */
