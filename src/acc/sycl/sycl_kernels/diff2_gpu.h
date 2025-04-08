#ifndef DIFF2_GPU_KERNELS_H_
#define DIFF2_GPU_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sycl/sycl.hpp>

#include "src/acc/acc_projector.h"
#include "src/acc/sycl/sycl_settings.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"
#include "src/acc/sycl/sycl_kernels/helper.h"
#include "src/acc/sycl/sycl_kernels/helper_gpu.h"

namespace syclGpuKernels
{

/*
 *   	DIFFERENCE-BASED KERNELS
 */
template<bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
void sycl_kernel_diff2_coarse(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_eulers, XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
		XFLOAT *g_real,	XFLOAT *g_imag, XFLOAT *g_corr, XFLOAT *g_diff2s,
		int trans_num, int image_size,
		XFLOAT *s_eulers, XFLOAT *s_ref_real, XFLOAT *s_ref_imag,
		XFLOAT *s_real, XFLOAT *s_imag, XFLOAT *s_corr
		)
{
	const int tid = nit.get_local_id(2);
	const int blockid = nit.get_group_linear_id();

	const int xSize = projector.imgX;
	const int ySize = projector.imgY;
	const int xySize = xSize * ySize;
	const int zSize = projector.imgZ;
	const int maxR = projector.maxR;

	const int max_block_pass_euler {((eulers_per_block*9)/block_sz + 1) * block_sz};
	for (int i = tid; i < max_block_pass_euler; i += block_sz)
		if (i < eulers_per_block * 9)
			s_eulers[i] = g_eulers[blockid*eulers_per_block*9 + i];

	XFLOAT diff2s[eulers_per_block] {0.0f};

	const XFLOAT tx {trans_x[tid % trans_num]};
	const XFLOAT ty {trans_y[tid % trans_num]};
	const XFLOAT tz {trans_z[tid % trans_num]};

	//Step through data
	const int max_block_pass_pixel {(image_size/block_sz + 1) * block_sz};
	__group_barrier(nit);
	for (int init_pixel = 0; init_pixel < max_block_pass_pixel; init_pixel += block_sz/prefetch_fraction)
	{

		//Prefetch block-fraction-wise
		if (init_pixel + tid/prefetch_fraction < image_size)
		{
			int x, y, z, xy;
			if (DATA3D)
			{
				z  = (init_pixel + tid/prefetch_fraction) / xySize;
				xy = (init_pixel + tid/prefetch_fraction) % xySize;
				x  = xy % xSize;
				y  = xy / xSize;
				if (z > maxR)
					z -= zSize;
			}
			else
			{
				x = (init_pixel + tid/prefetch_fraction) % xSize;
				y = (init_pixel + tid/prefetch_fraction) / xSize;
			}
			if (y > maxR)
				y -= ySize;

			for (int e = tid%prefetch_fraction; e < eulers_per_block; e += prefetch_fraction)
			{
				if (DATA3D) // if DATA3D, then REF3D as well.
				{
					projector.project3Dmodel(
						x, y, z,
						s_eulers[e*9  ], s_eulers[e*9+1], s_eulers[e*9+2],
						s_eulers[e*9+3], s_eulers[e*9+4], s_eulers[e*9+5],
						s_eulers[e*9+6], s_eulers[e*9+7], s_eulers[e*9+8],
						s_ref_real[eulers_per_block * (tid/prefetch_fraction) + e],
						s_ref_imag[eulers_per_block * (tid/prefetch_fraction) + e]);
				}
				else if (REF3D)
				{
					projector.project3Dmodel(
						x,y,
						s_eulers[e*9  ], s_eulers[e*9+1], s_eulers[e*9+3],
						s_eulers[e*9+4], s_eulers[e*9+6], s_eulers[e*9+7],
						s_ref_real[eulers_per_block * (tid/prefetch_fraction) + e],
						s_ref_imag[eulers_per_block * (tid/prefetch_fraction) + e]);
				}
				else
				{
					projector.project2Dmodel(
						x,y,
						s_eulers[e*9  ], s_eulers[e*9+1], s_eulers[e*9+3], s_eulers[e*9+4],
						s_ref_real[eulers_per_block * (tid/prefetch_fraction) + e],
						s_ref_imag[eulers_per_block * (tid/prefetch_fraction) + e]);
				}
			}
		}

		//Prefetch block-wise
		if (init_pixel % block_sz == 0 && init_pixel + tid < image_size)
		{
			s_real[tid] = g_real[init_pixel + tid];
			s_imag[tid] = g_imag[init_pixel + tid];
			s_corr[tid] = g_corr[init_pixel + tid] / 2.0f;
		}
		__group_barrier(nit);

		if (tid/trans_num < block_sz/trans_num) // NOTE int division A/B==C/B !=> A==C
		{
			for (int pix = tid/trans_num; pix < block_sz/prefetch_fraction; pix += block_sz/trans_num)
			{
				if ((init_pixel + pix) >= image_size) break;

				int x, y, z, xy;
				if (DATA3D)
				{
					z  = (init_pixel + pix) / xySize;
					xy = (init_pixel + pix) % xySize;
					x  = xy % xSize;
					y  = xy / ySize;
					if (z > maxR)
						z -= zSize;
				}
				else
				{
					x = (init_pixel + pix) % xSize;
					y = (init_pixel + pix) / xSize;
				}
				if (y > maxR)
					y -= ySize;

				XFLOAT real, imag;
				if (DATA3D)
				{
//					translatePixel(x, y, z, tx, ty, tz, s_real[pix+init_pixel%block_sz], s_imag[pix+init_pixel%block_sz], real, imag);
					XFLOAT val {x*tx + y*ty + z*tz};
#ifdef ACC_DOUBLE_PRECISION
					XFLOAT s {sycl::sin(val)};
					XFLOAT c {sycl::cos(val)};
#else
					XFLOAT s {sycl::native::sin(val)};
					XFLOAT c {sycl::native::cos(val)};
#endif
					real = c * s_real[pix + init_pixel%block_sz] - s * s_imag[pix + init_pixel%block_sz];
					imag = c * s_imag[pix + init_pixel%block_sz] + s * s_real[pix + init_pixel%block_sz];
				}
				else
				{
//					translatePixel(x, y,    tx, ty,     s_real[pix+init_pixel%block_sz], s_imag[pix+init_pixel%block_sz], real, imag);
					XFLOAT val {x*tx + y*ty};
#ifdef ACC_DOUBLE_PRECISION
					XFLOAT s {sycl::sin(val)};
					XFLOAT c {sycl::cos(val)};
#else
					XFLOAT s {sycl::native::sin(val)};
					XFLOAT c {sycl::native::cos(val)};
#endif
					real = c*s_real[pix + init_pixel%block_sz] - s*s_imag[pix + init_pixel%block_sz];
					imag = c*s_imag[pix + init_pixel%block_sz] + s*s_real[pix + init_pixel%block_sz];
				}

				for (int e = 0; e < eulers_per_block; e++)
				{
					XFLOAT diff_real {s_ref_real[eulers_per_block*pix + e] - real};
					XFLOAT diff_imag {s_ref_imag[eulers_per_block*pix + e] - imag};
					diff2s[e] += (diff_real*diff_real + diff_imag*diff_imag) * s_corr[pix + init_pixel%block_sz];
				}
			}	// for over i = tid/trans_num
		}	// if < image_size
	}	// for over init_pixel

	//Set global
	for (int e = 0; e < eulers_per_block; e++)
		atomic_ref_wg( g_diff2s[blockid*eulers_per_block*trans_num + e*trans_num + tid%trans_num] ).fetch_add(diff2s[e]);
}

template <bool REF3D, bool DATA3D, int block_sz>
void sycl_kernel_diff2_fine(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		int image_size, XFLOAT sum_init, int todo_blocks,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		XFLOAT *s)
{
	const int bid = nit.get_group_linear_id();
	const int tid = nit.get_local_id(2);

	const int xSize = projector.imgX;
	const int ySize = projector.imgY;
	const int zSize = projector.imgZ;
	const int maxR = projector.maxR;

	if (bid < todo_blocks ) // we only need to make
	{
		int trans_num = static_cast<int>(d_job_num[bid]); //how many transes we have for this rot
		for (int itrans = 0; itrans < trans_num; itrans++)
			s[itrans*block_sz + tid] = 0.0f;

		// index of comparison
		int ix = static_cast<int>(d_rot_idx[d_job_idx[bid]]);
		int pass_num {image_size/block_sz + 1};

		for (int pass = 0; pass < pass_num; pass++) // finish an entire ref image each block
		{
			int pixel = pass*block_sz + tid;

			if (pixel < image_size)
			{
				int x, y, z, xy;
				if (DATA3D)
				{
					z  = pixel / (xSize*ySize);
					xy = pixel % (xSize*ySize);
					x  = xy % xSize;
					y  = xy / xSize;
					if (z > maxR)
					{
						if (z >= zSize - maxR)
							z = z - zSize;
						else
							x = maxR;
					}
				}
				else
				{
					x = pixel % xSize;
					y = pixel / xSize;
				}

				if (y > maxR)
				{
					if (y >= ySize - maxR)
						y = y - ySize;
					else
						x = maxR;
				}

				XFLOAT ref_real, ref_imag;
				if (DATA3D)
					projector.project3Dmodel(
							x, y, z,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1], g_eulers[ix*9 + 2],
							g_eulers[ix*9 + 3], g_eulers[ix*9 + 4], g_eulers[ix*9 + 5],
							g_eulers[ix*9 + 6], g_eulers[ix*9 + 7], g_eulers[ix*9 + 8],
							ref_real, ref_imag);
				else if (REF3D)
					projector.project3Dmodel(
							x, y,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1], g_eulers[ix*9 + 3],
							g_eulers[ix*9 + 4], g_eulers[ix*9 + 6], g_eulers[ix*9 + 7],
							ref_real, ref_imag);
				else
					projector.project2Dmodel(
							x, y,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1],
							g_eulers[ix*9 + 3], g_eulers[ix*9 + 4],
							ref_real, ref_imag);

				for (int itrans = 0; itrans < trans_num; itrans++) // finish all translations in each partial pass
				{
					int iy = static_cast<int>(d_trans_idx[d_job_idx[bid]]) + itrans;

					XFLOAT shifted_real, shifted_imag;
					if (DATA3D)
						translatePixel(x, y, z, trans_x[iy], trans_y[iy], trans_z[iy], g_imgs_real[pixel], g_imgs_imag[pixel], shifted_real, shifted_imag);
					else
						translatePixel(x, y,    trans_x[iy], trans_y[iy],              g_imgs_real[pixel], g_imgs_imag[pixel], shifted_real, shifted_imag);

					XFLOAT diff_real =  ref_real - shifted_real;
					XFLOAT diff_imag =  ref_imag - shifted_imag;
					s[itrans*block_sz + tid] += (diff_real*diff_real + diff_imag*diff_imag) * 0.5f * g_corr_img[pixel];
				}
			}
		}
		__group_barrier(nit);

// TODO: reduction can be used
		for (int j = block_sz/2; j > 0; j /= 2)
		{
			if (tid < j)
				for (int itrans = 0; itrans < trans_num; itrans++) // finish all translations in each partial pass
					s[itrans*block_sz + tid] += s[itrans*block_sz + tid + j];

			__group_barrier(nit);
		}

		if (tid < trans_num)
			g_diff2s[d_job_idx[bid] + tid] += s[tid * block_sz] + sum_init;
	}
}


/*
 *   	CROSS-CORRELATION-BASED KERNELS
 */
template<bool REF3D, bool DATA3D, int block_sz>
void sycl_kernel_diff2_CC_coarse(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_corr_img, XFLOAT *g_diff2,
		int trans_num, int image_size,
		XFLOAT *s_weight, XFLOAT *s_norm
		)
{
	const int iorient = nit.get_group(0);
	const int itrans = nit.get_group(1);
	const int tid = nit.get_local_id(2);

	const int xSize = projector.imgX;
	const int ySize = projector.imgY;
	const int zSize = projector.imgZ;
	const int maxR = projector.maxR;

	const XFLOAT e0 = g_eulers[iorient * 9    ];
	const XFLOAT e1 = g_eulers[iorient * 9 + 1];
	const XFLOAT e2 = g_eulers[iorient * 9 + 2];
	const XFLOAT e3 = g_eulers[iorient * 9 + 3];
	const XFLOAT e4 = g_eulers[iorient * 9 + 4];
	const XFLOAT e5 = g_eulers[iorient * 9 + 5];
	const XFLOAT e6 = g_eulers[iorient * 9 + 6];
	const XFLOAT e7 = g_eulers[iorient * 9 + 7];
	const XFLOAT e8 = g_eulers[iorient * 9 + 8];

	s_weight[tid] = 0.0f;
	s_norm[tid] = 0.0f;

	const int pixel_pass_num {image_size/block_sz + 1};
	for (int pass = 0; pass < pixel_pass_num; pass++)
	{
		const int pixel {pass*block_sz + tid};
		if (pixel < image_size)
		{
			int x, y, z, xy;
			if (DATA3D)
			{
				z =  pixel / (xSize*ySize);
				xy = pixel % (xSize*ySize);
				x =  xy % xSize;
				y =  xy / xSize;
				if (z > maxR)
				{
					if (z >= zSize - maxR)
						z = z - zSize;
					else
						x = maxR;
				}
			}
			else
			{
				x = pixel % xSize;
				y = pixel / xSize;
			}
			if (y > maxR)
			{
				if (y >= ySize - maxR)
					y = y - ySize;
				else
					x = maxR;
			}

			XFLOAT ref_real, ref_imag;
			if (DATA3D)
				projector.project3Dmodel(x, y, z, e0, e1, e2, e3, e4, e5, e6, e7, e8, ref_real, ref_imag);
			else if (REF3D)
				projector.project3Dmodel(x, y, e0, e1, e3, e4, e6, e7, ref_real, ref_imag);
			else
				projector.project2Dmodel(x, y, e0, e1, e3, e4, ref_real, ref_imag);

			XFLOAT real, imag;
			if (DATA3D)
				translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], g_imgs_real[pixel], g_imgs_imag[pixel], real, imag);
			else
				translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    g_imgs_real[pixel], g_imgs_imag[pixel], real, imag);

			s_weight[tid] += (ref_real *     real + ref_imag *     imag) * g_corr_img[pixel];
			s_norm  [tid] += (ref_real * ref_real + ref_imag * ref_imag) * g_corr_img[pixel];
		}
	}
	__group_barrier(nit);

// TODO: This could be optimized using sub_group collective and s_weight and s_norm can be thread private
	for (int j = block_sz/2; j > 0; j /= 2)
	{
		if (tid < j)
		{
			s_weight[tid] += s_weight[tid+j];
			s_norm[tid]   += s_norm[tid+j];
		}
		__group_barrier(nit);
	}

	if (tid == 0)
#ifdef ACC_DOUBLE_PRECISION
		g_diff2[iorient*trans_num + itrans] += -s_weight[0] / sycl::sqrt(s_norm[0]);
#else
		g_diff2[iorient*trans_num + itrans] += -s_weight[0] / sycl::native::sqrt(s_norm[0]);
#endif
}

template<bool REF3D, bool DATA3D, int block_sz>
void sycl_kernel_diff2_CC_fine(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		int image_size, XFLOAT sum_init, int todo_blocks,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		XFLOAT *s, XFLOAT *s_cc)
{
	const int bid = nit.get_group_linear_id();
	const int tid = nit.get_local_id(2);

	const int xSize = projector.imgX;
	const int ySize = projector.imgY;
	const int zSize = projector.imgZ;
	const int maxR = projector.maxR;

	if (bid < todo_blocks) // we only need to make
	{
		int trans_num = static_cast<int>(d_job_num[bid]); //how many transes we have for this rot
		for (int itrans = 0; itrans < trans_num; itrans++)
		{
			s[   itrans*block_sz + tid] = 0.0f;
			s_cc[itrans*block_sz + tid] = 0.0f;
		}

		// index of comparison
		int ix = static_cast<int>(d_rot_idx[d_job_idx[bid]]);
		int pass_num {image_size/block_sz + 1};
		for (int pass = 0; pass < pass_num; pass++) // finish an entire ref image each block
		{
			int pixel = pass*block_sz + tid;

			if (pixel < image_size)
			{
				int x, y, z, xy;
				if (DATA3D)
				{
					z  = pixel / (xSize*ySize);
					xy = pixel % (xSize*ySize);
					x  = xy % xSize;
					y  = xy / xSize;
					if (z > maxR)
					{
						if (z >= zSize - maxR)
							z = z - zSize;
						else
							x = maxR;
					}
				}
				else
				{
					x = pixel % xSize;
					y = pixel / xSize;
				}

				if (y > maxR)
				{
					if (y >= ySize - maxR)
						y = y - ySize;
					else
						x = maxR;
				}

				XFLOAT ref_real, ref_imag;
				if (DATA3D)
					projector.project3Dmodel(
							x, y, z,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1], g_eulers[ix*9 + 2],
							g_eulers[ix*9 + 3], g_eulers[ix*9 + 4], g_eulers[ix*9 + 5],
							g_eulers[ix*9 + 6], g_eulers[ix*9 + 7], g_eulers[ix*9 + 8],
							ref_real, ref_imag);
				else if (REF3D)
					projector.project3Dmodel(
							x, y,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1], g_eulers[ix*9 + 3],
							g_eulers[ix*9 + 4], g_eulers[ix*9 + 6], g_eulers[ix*9 + 7],
							ref_real, ref_imag);
				else
					projector.project2Dmodel(
							x, y,
							g_eulers[ix*9    ], g_eulers[ix*9 + 1],
							g_eulers[ix*9 + 3], g_eulers[ix*9 + 4],
							ref_real, ref_imag);

				for (int itrans = 0; itrans < trans_num; itrans++) // finish all translations in each partial pass
				{
					int iy = static_cast<int>(d_trans_idx[d_job_idx[bid]]) + itrans;

					XFLOAT shifted_real, shifted_imag;
					if(DATA3D)
						translatePixel(x, y, z, g_trans_x[iy], g_trans_y[iy], g_trans_z[iy], g_imgs_real[pixel], g_imgs_imag[pixel], shifted_real, shifted_imag);
					else
						translatePixel(x, y,    g_trans_x[iy], g_trans_y[iy],                g_imgs_real[pixel], g_imgs_imag[pixel], shifted_real, shifted_imag);

					s   [itrans*block_sz + tid] += (ref_real*shifted_real + ref_imag*shifted_imag) * g_corr_img[pixel];
					s_cc[itrans*block_sz + tid] += (ref_real*    ref_real + ref_imag*    ref_imag) * g_corr_img[pixel];
				}
			}
		}
		__group_barrier(nit);

		for (int j = block_sz/2; j > 0; j /= 2)
		{
			if (tid < j)
			{
				for (int itrans = 0; itrans < trans_num; itrans++) // finish all translations in each partial pass
				{
					s[   itrans*block_sz + tid] += s[   itrans*block_sz + tid + j];
					s_cc[itrans*block_sz + tid] += s_cc[itrans*block_sz + tid + j];
				}
			}
			__group_barrier(nit);
		}

		if (tid < trans_num)
#ifdef ACC_DOUBLE_PRECISION
			g_diff2s[d_job_idx[bid] + tid] += -s[tid * block_sz] / sycl::sqrt(s_cc[tid * block_sz]);
#else
			g_diff2s[d_job_idx[bid] + tid] += -s[tid * block_sz] / sycl::native::sqrt(s_cc[tid * block_sz]);
#endif
    }
}

} // end of namespace syclKernels

#endif /* DIFF2_GPU_KERNELS_H_ */
