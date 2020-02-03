#ifndef WAVG_KERNEL_H_
#define WAVG_KERNEL_H_

#include <vector>
#include <iostream>
#include <fstream>

#include "src/acc/cpu/cpu_settings.h"
#include "src/acc/acc_projector.h"
#include "src/acc/cpu/cpu_kernels/cpu_utils.h"
#include "src/acc/cpu/cpu_kernels/helper.h"

namespace CpuKernels
{

// sincos lookup table optimization. Function translatePixel calls
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions.
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty).
template<bool CTFPREMULTIPLIED, bool REFCTF, bool REF3D>
void wavg_ref3D(
		XFLOAT * RESTRICT   g_eulers,
		AccProjectorKernel &projector,
		unsigned long       image_size,
		unsigned long       orientation_num,
#ifdef DEBUG_CUDA
		XFLOAT * RESTRICT   _g_img_real,
#else
		XFLOAT * RESTRICT   g_img_real,
#endif
		XFLOAT * RESTRICT   g_img_imag,
		XFLOAT * RESTRICT   g_trans_x,
		XFLOAT * RESTRICT   g_trans_y,
		XFLOAT * RESTRICT   g_trans_z,
		XFLOAT * RESTRICT   g_weights,
		XFLOAT * RESTRICT   g_ctfs,
		XFLOAT * RESTRICT   g_wdiff2s_parts,
		XFLOAT * RESTRICT   g_wdiff2s_AA,
		XFLOAT * RESTRICT   g_wdiff2s_XA,
		unsigned long       trans_num,
		XFLOAT              weight_norm,
		XFLOAT              significant_weight,
		XFLOAT              part_scale)
{
#ifdef DEBUG_CUDA
	checkedArray<XFLOAT> g_img_real;
	g_img_real.initCheckedArray(_g_img_real);
#endif
	XFLOAT ref_real, ref_imag, img_real, img_imag, trans_real, trans_imag;

	for(unsigned long bid=0; bid<orientation_num; bid++) {

		// Copy the rotation matrix to local variables
		int offset = bid * 9;
		XFLOAT e0 = g_eulers[offset  ], e1 = g_eulers[offset+1];
		XFLOAT e3 = g_eulers[offset+3], e4 = g_eulers[offset+4];
		XFLOAT e6 = g_eulers[offset+6], e7 = g_eulers[offset+7];

		// pre-compute sin and cos for x and y direction
		int xSize = projector.imgX;
		int ySize = projector.imgY;
		XFLOAT sin_x[trans_num][xSize], cos_x[trans_num][xSize];
		XFLOAT sin_y[trans_num][ySize], cos_y[trans_num][ySize];

		computeSincosLookupTable2D(trans_num, g_trans_x, g_trans_y, xSize, ySize,
								  &sin_x[0][0], &cos_x[0][0],
								  &sin_y[0][0], &cos_y[0][0]);

		XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
		unsigned long pixel = 0;
		for(int iy = 0; iy < ySize; iy++) {
			int xstart = 0, xend = xSize;
			int y = iy;
			if (iy > projector.maxR) {
				if (iy >= ySize - projector.maxR)
					y = iy - ySize;
				else {
					// handle special case for one pixel
					xstart = projector.maxR;
					xend   = xstart + 1;
				}
			}

			XFLOAT ref_real[xSize], ref_imag[xSize];
			XFLOAT img_real[xSize], img_imag[xSize];

			#pragma omp simd
			for(int x = xstart; x < xend; x++) {
				if(REF3D)
					projector.project3Dmodel(x, y, e0, e1, e3, e4, e6, e7,
											ref_real[x], ref_imag[x]);
				else
					projector.project2Dmodel(x, y, e0, e1, e3, e4,
											ref_real[x], ref_imag[x]);
				if (REFCTF)
				{
					if(CTFPREMULTIPLIED)
					{
						ref_real[x] *= g_ctfs[pixel + x] * g_ctfs[pixel + x];
						ref_imag[x] *= g_ctfs[pixel + x] * g_ctfs[pixel + x];
					}
					else
					{
						ref_real[x] *= g_ctfs[pixel + x];
						ref_imag[x] *= g_ctfs[pixel + x];
					}
				}
				else {
					ref_real[x] *= part_scale;
					ref_imag[x] *= part_scale;
				}

				img_real[x] = g_img_real[pixel + x];
				img_imag[x] = g_img_imag[pixel + x];
			}

			for (unsigned long itrans = 0; itrans < trans_num; itrans++) {
				XFLOAT weight = g_weights[bid * trans_num + itrans];
				if (weight < significant_weight)
					continue;

				weight *= weight_norm_inverse;
				XFLOAT trans_cos_y, trans_sin_y;
				if ( y < 0) {
					trans_cos_y =  cos_y[itrans][-y];
					trans_sin_y = -sin_y[itrans][-y];
				}
				else {
					trans_cos_y = cos_y[itrans][y];
					trans_sin_y = sin_y[itrans][y];
				}

				XFLOAT *trans_cos_x = &cos_x[itrans][0];
				XFLOAT *trans_sin_x = &sin_x[itrans][0];

				for(int x = xstart; x < xend; x++) {

					XFLOAT ss = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
					XFLOAT cc = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

					XFLOAT trans_real = cc * img_real[x] - ss * img_imag[x];
					XFLOAT trans_imag = cc * img_imag[x] + ss * img_real[x];
					/*
					XFLOAT trans_real, trans_imag;
					translatePixel(x, y,  g_trans_x[itrans], g_trans_y[itrans],
								   img_real[x], img_imag[x], trans_real, trans_imag);
					*/
					XFLOAT diff_real = ref_real[x] - trans_real;
					XFLOAT diff_imag = ref_imag[x] - trans_imag;

					g_wdiff2s_parts[pixel + x] += weight * (diff_real    * diff_real   + diff_imag    * diff_imag);
					g_wdiff2s_XA   [pixel + x] += weight * (ref_real[x]  * trans_real  + ref_imag[x]  * trans_imag);
					g_wdiff2s_AA   [pixel + x] += weight * (ref_real[x]  * ref_real[x] + ref_imag[x]  * ref_imag[x] );
				}
			}  // for itrans

			pixel += (unsigned long)xSize;
		} // y direction
	} // bid
}

template<bool CTFPREMULTIPLIED, bool REFCTF>
void wavg_3D(
		XFLOAT * RESTRICT   g_eulers,
		AccProjectorKernel &projector,
		unsigned long       image_size,
		unsigned long       orientation_num,
#ifdef DEBUG_CUDA
		XFLOAT * RESTRICT   _g_img_real,
#else
		XFLOAT * RESTRICT   g_img_real,
#endif
		XFLOAT * RESTRICT   g_img_imag,
		XFLOAT * RESTRICT   g_trans_x,
		XFLOAT * RESTRICT   g_trans_y,
		XFLOAT * RESTRICT   g_trans_z,
		XFLOAT * RESTRICT   g_weights,
		XFLOAT * RESTRICT   g_ctfs,
		XFLOAT * RESTRICT   g_wdiff2s_parts,
		XFLOAT * RESTRICT   g_wdiff2s_AA,
		XFLOAT * RESTRICT   g_wdiff2s_XA,
		unsigned long       trans_num,
		XFLOAT              weight_norm,
		XFLOAT              significant_weight,
		XFLOAT              part_scale)
{
#ifdef DEBUG_CUDA
	checkedArray<XFLOAT> g_img_real;
	g_img_real.initCheckedArray(_g_img_real);
#endif
	XFLOAT ref_real, ref_imag, img_real, img_imag, trans_real, trans_imag;

	for(unsigned long bid=0; bid<orientation_num; bid++) {
		// Copy the rotation matrix to local variables
		int offset = bid * 9;
		XFLOAT e0 = g_eulers[offset  ], e1 = g_eulers[offset+1];
		XFLOAT e2 = g_eulers[offset+2], e3 = g_eulers[offset+3];
		XFLOAT e4 = g_eulers[offset+4], e5 = g_eulers[offset+5];
		XFLOAT e6 = g_eulers[offset+6], e7 = g_eulers[offset+7];
		XFLOAT e8 = g_eulers[offset+8];

		// pre-compute sin and cos for x and y direction
		int xSize = projector.imgX;
		int ySize = projector.imgY;
		int zSize = projector.imgZ;
		XFLOAT sin_x[trans_num][xSize], cos_x[trans_num][xSize];
		XFLOAT sin_y[trans_num][ySize], cos_y[trans_num][ySize];
		XFLOAT sin_z[trans_num][zSize], cos_z[trans_num][zSize];

		computeSincosLookupTable3D(trans_num, g_trans_x, g_trans_y, g_trans_z,
								   xSize, ySize, zSize,
								  &sin_x[0][0], &cos_x[0][0],
								  &sin_y[0][0], &cos_y[0][0],
								  &sin_z[0][0], &cos_z[0][0]);

		XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
		unsigned long pixel = 0;
		for(int iz = 0; iz < zSize; iz ++) {
			int xstart_z = 0, xend_z = xSize;
			int z = iz;
			if (z > projector.maxR)
			{
				if (z >= zSize - projector.maxR)
					z = z - projector.imgZ;
				else
					xstart_z = projector.maxR;
					xend_z   = xstart_z + 1;
			}

			for(int iy = 0; iy < ySize; iy++) {
				int xstart_y = xstart_z, xend_y = xend_z;
				int y = iy;
				if (iy > projector.maxR) {
					if (iy >= ySize - projector.maxR)
						y = iy - ySize;
					else {
						xstart_y = projector.maxR;
						xend_y   = xstart_y + 1;
					}
				}

				XFLOAT ref_real[xSize],  ref_imag[xSize];
				XFLOAT img_real[xSize], img_imag[xSize];

				#pragma omp simd
				for(int x = xstart_y; x < xend_y; x++) {
					projector.project3Dmodel(x, y, z, e0, e1, e2, e3, e4, e5, e6, e7, e8,
											 ref_real[x], ref_imag[x]);
					if (REFCTF)
					{
						if(CTFPREMULTIPLIED)
						{
							ref_real[x] *= g_ctfs[pixel + x] * g_ctfs[pixel + x];
							ref_imag[x] *= g_ctfs[pixel + x] * g_ctfs[pixel + x];
						}
						else
						{
							ref_real[x] *= g_ctfs[pixel + x];
							ref_imag[x] *= g_ctfs[pixel + x];
						}
					}
					else {
						ref_real[x] *= part_scale;
						ref_imag[x] *= part_scale;
					}

					img_real[x] = g_img_real[pixel + x];
					img_imag[x] = g_img_imag[pixel + x];
				}

				for (unsigned long itrans = 0; itrans < trans_num; itrans++) {
					XFLOAT weight = g_weights[bid * trans_num + itrans];
					if (weight < significant_weight)
						continue;

					weight *= weight_norm_inverse;
					XFLOAT trans_cos_z, trans_sin_z;
					if ( z < 0) {
						trans_cos_z =  cos_z[itrans][-z];
						trans_sin_z = -sin_z[itrans][-z];
					}
					else {
						trans_cos_z = cos_z[itrans][z];
						trans_sin_z = sin_z[itrans][z];
					}

					XFLOAT trans_cos_y, trans_sin_y;
					if ( y < 0) {
						trans_cos_y =  cos_y[itrans][-y];
						trans_sin_y = -sin_y[itrans][-y];
					}
					else {
						trans_cos_y = cos_y[itrans][y];
						trans_sin_y = sin_y[itrans][y];
					}

					XFLOAT *trans_cos_x = &cos_x[itrans][0];
					XFLOAT *trans_sin_x = &sin_x[itrans][0];

					for(int x = xstart_y; x < xend_y; x++) {
						// TODO check the math
						XFLOAT s  = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
						XFLOAT c  = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

						XFLOAT ss = s * trans_cos_z + c * trans_sin_z;
						XFLOAT cc = c * trans_cos_z - s * trans_sin_z;

						XFLOAT trans_real = cc * img_real[x] - ss * img_imag[x];
						XFLOAT trans_imag = cc * img_imag[x] + ss * img_real[x];
						/*
						XFLOAT trans_real, trans_imag;
						translatePixel(x, y,  g_trans_x[itrans], g_trans_y[itrans],
									   img_real[x], img_imag[x], trans_real, trans_imag);

						where translatePixel is:
							sincosf( x * tx + y * ty , &s, &c );
							tReal = c * real - s * imag;
							tImag = c * imag + s * real;
						*/
						XFLOAT diff_real = ref_real[x] - trans_real;
						XFLOAT diff_imag = ref_imag[x] - trans_imag;

						g_wdiff2s_parts[pixel + x] += weight * (diff_real    * diff_real   + diff_imag    * diff_imag);
						g_wdiff2s_XA   [pixel + x] += weight * (ref_real[x]  * trans_real  + ref_imag[x]  * trans_imag);
						g_wdiff2s_AA   [pixel + x] += weight * (ref_real[x]  * ref_real[x] + ref_imag[x]  * ref_imag[x] );
					}
				}  // for itrans

				pixel += (unsigned long)xSize;
			} // y direction
		} // z direction
	} // bid
}

} // end of namespace CpuKernels

#endif /* WAVG_KERNEL_H_ */
