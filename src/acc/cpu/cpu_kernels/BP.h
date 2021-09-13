#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <cassert>

#include "src/acc/acc_backprojector.h"
#include "src/acc/cpu/cpu_kernels/helper.h"

namespace CpuKernels
{
template < bool CTF_PREMULTIPLIED >
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject2D(
		unsigned long imageCount,
		int     block_size,
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
		int mdl_inity,
		tbb::spin_mutex *mutexes)
{
	int img_y_half = img_y / 2;

	int max_r2_out = max_r2 * padding_factor * padding_factor;

	// pre-compute sin and cos for x and y direction
	XFLOAT sin_x[translation_num][img_x], cos_x[translation_num][img_x];
	XFLOAT sin_y[translation_num][img_y], cos_y[translation_num][img_y];

	computeSincosLookupTable2D(translation_num, g_trans_x, g_trans_y, 
								img_x, img_y,
								&sin_x[0][0], &cos_x[0][0], 
								&sin_y[0][0], &cos_y[0][0]);

	// Set up some other variables
	XFLOAT s_eulers[4];
	
	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
	
	XFLOAT xp[img_x], yp[img_x];
	XFLOAT real[img_x], imag[img_x], Fweight[img_x];
	
	for (unsigned long img=0; img<imageCount; img++) {

		// Copy the rotation matrix to local variables
		s_eulers[0] = g_eulers[img*9+0] * padding_factor;
		s_eulers[1] = g_eulers[img*9+1] * padding_factor;
		s_eulers[2] = g_eulers[img*9+3] * padding_factor;
		s_eulers[3] = g_eulers[img*9+4] * padding_factor;
		
		size_t pixel = 0;
		
		for(int iy = 0; iy < img_y; iy++) {
			int y = iy;
			if (iy > img_y_half) {
				y = iy - img_y;
			}
			
			int xmax = img_x;

			memset(Fweight,0,sizeof(XFLOAT)*img_x);
			memset(real,   0,sizeof(XFLOAT)*img_x);
			memset(imag,   0,sizeof(XFLOAT)*img_x);
			
			for (unsigned long itrans = 0; itrans < translation_num; itrans++)
			{
				XFLOAT weight = g_weights[img * translation_num + itrans];
				
				if (weight < significant_weight) 
					continue;

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
			
				for(int x=0; x<xmax; x++) {
					//WAVG
					XFLOAT minvsigma2 = g_Minvsigma2s[pixel + x];
					XFLOAT ctf = g_ctfs[pixel + x];
					XFLOAT my_weight;
					
					if(CTF_PREMULTIPLIED) {
						my_weight = weight * weight_norm_inverse * minvsigma2;
						Fweight[x] += my_weight  * ctf;
					}
					else {
						my_weight = weight * weight_norm_inverse * ctf * minvsigma2;
						Fweight[x] += my_weight  * ctf;
					}
					/*
					CpuKernels::translatePixel(x, y, 
					 * g_trans_x[itrans], g_trans_y[itrans], 
					 * img_real, img_imag, temp_real, temp_imag);
				     */
					XFLOAT img_real = g_img_real[pixel + x];
					XFLOAT img_imag = g_img_imag[pixel + x];
				
					XFLOAT ss = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
					XFLOAT cc = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

					XFLOAT temp_real = cc * img_real - ss * img_imag;
					XFLOAT temp_imag = cc * img_imag + ss * img_real;

					real[x] += temp_real * my_weight;
					imag[x] += temp_imag * my_weight;
				}  // for x
			}  // for itrans

			for(int x=0; x<xmax; x++) {	
				if (Fweight[x] <= (XFLOAT) 0.0)
					continue;
				
				// Get logical coordinates in the 3D map
				xp[x] = (s_eulers[0] * x + s_eulers[1] * y );
				yp[x] = (s_eulers[2] * x + s_eulers[3] * y );
				
				// Only consider pixels that are projected inside the allowed circle in output coordinates.
				//     --JZ, Nov. 26th 2018
				if ( ( xp[x] * xp[x] + yp[x] * yp[x] ) > max_r2_out)
				{
					Fweight[x]= (XFLOAT) 0.0;
					continue;
				}

				// Only asymmetric half is stored
				if (xp[x] < (XFLOAT) 0.0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp[x] = -xp[x];
					yp[x] = -yp[x];
					imag[x] = -imag[x];
				}
			}  // for x
			
			for(int x=0; x<xmax; x++) {
				if (Fweight[x] <= (XFLOAT) 0.0)
					continue;
				
				int x0 = floorf(xp[x]);
				XFLOAT fx = xp[x] - x0;

				int y0 = floorf(yp[x]);
				XFLOAT fy = yp[x] - y0;
				y0 -= mdl_inity;
				int y1 = y0 + 1;
				
				XFLOAT mfx = (XFLOAT) 1.0 - fx;
				XFLOAT mfy = (XFLOAT) 1.0 - fy;

				XFLOAT dd00 = mfy * mfx;
				XFLOAT dd01 = mfy *  fx;
				XFLOAT dd10 =  fy * mfx;
				XFLOAT dd11 =  fy *  fx;

				size_t idx_tmp;

				// Locking necessary since all threads share the same back projector
				{
					idx_tmp = (size_t)y0 * (size_t)mdl_x + (size_t)x0;

					tbb::spin_mutex::scoped_lock lock(mutexes[y0]);
					g_model_real  [idx_tmp] += dd00 * real[x];
					g_model_imag  [idx_tmp] += dd00 * imag[x];
					g_model_weight[idx_tmp] += dd00 * Fweight[x];

					idx_tmp = idx_tmp + 1;  // x1 = x0 + 1

					g_model_real  [idx_tmp] += dd01 * real[x];
					g_model_imag  [idx_tmp] += dd01 * imag[x];
					g_model_weight[idx_tmp] += dd01 * Fweight[x];
				}  // scoping for first lock

				{
					idx_tmp = (size_t)y1 * (size_t)mdl_x + (size_t)x0;

					tbb::spin_mutex::scoped_lock lock(mutexes[y1]);
					g_model_real  [idx_tmp] += dd10 * real[x];
					g_model_imag  [idx_tmp] += dd10 * imag[x];
					g_model_weight[idx_tmp] += dd10 * Fweight[x];

					idx_tmp = idx_tmp + 1;  // x1 = x0 + 1
					g_model_real  [idx_tmp] += dd11 * real[x];
					g_model_imag  [idx_tmp] += dd11 * imag[x];
					g_model_weight[idx_tmp] += dd11 * Fweight[x];
				}  // scoping for second lock
			}  // for x
			
			pixel += (size_t)img_x;
		} // for y
	} // for img
}

template < bool DATA3D, bool CTF_PREMULTIPLIED >
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject3D(
		unsigned long imageCount,
		int     block_size,
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
		size_t   img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz,
		tbb::spin_mutex *mutexes)
{
	int img_y_half = img_y / 2;
	int img_z_half = img_z / 2;

	int max_r2_vol = max_r2 * padding_factor * padding_factor;
	
	// Set up some variables
	XFLOAT s_eulers[9];
	
	//   We collect block_size number of values before storing the results to
	//   help vectorization and control memory accesses
	XFLOAT real[block_size], imag[block_size], Fweight[block_size];
	XFLOAT xp[block_size], yp[block_size], zp[block_size];

	for (unsigned long img=0; img<imageCount; img++) {

		 for (int i = 0; i < 9; i++)
			 s_eulers[i] = g_eulers[img*9+i];

		XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;

		int pixel_pass_num(0);
		if(DATA3D)
			pixel_pass_num = (ceilf((float)img_xyz/(float)block_size));
		else
			pixel_pass_num = (ceilf((float)img_xyz/(float)block_size));

		for (unsigned pass = 0; pass < pixel_pass_num; pass++)
		{
			memset(Fweight,0,sizeof(XFLOAT)*block_size);
			#pragma omp simd
			for(int tid=0; tid<block_size; tid++)
			{
				int ok_for_next(1);  // This flag avoids continues, helping the vectorizer

				size_t pixel(0);
				if(DATA3D)
					pixel = ((size_t)pass * (size_t)block_size) + (size_t)tid;
				else
					pixel = ((size_t)pass * (size_t)block_size) + (size_t)tid;

				if (pixel >= img_xyz)
					continue; // just doesn't make sense to proceed in this case

				int x,y,z,xy;
				XFLOAT minvsigma2, ctf, img_real, img_imag, weight;

				if(DATA3D)
				{
					z =  CpuKernels::floorfracf(pixel, (size_t)img_x*(size_t)img_y);
					xy = pixel % (img_x*img_y);
					x =             xy  % img_x;
					y = CpuKernels::floorfracf( xy,   (size_t)img_x);

					if (z > img_z_half)
					{
						z = z - img_z;

						if(x==0)
							ok_for_next=0;
					}
				}
				else
				{
					x =             pixel % img_x;
					y = CpuKernels::floorfracf( pixel , (size_t)img_x);
				}
				if (y > img_y_half)
				{
					y = y - img_y;
				}
				// Get logical coordinates in the 3D map

				if(DATA3D)
				{
					xp[tid] = (s_eulers[0] * x + s_eulers[1] * y + s_eulers[2] * z) * padding_factor;
					yp[tid] = (s_eulers[3] * x + s_eulers[4] * y + s_eulers[5] * z) * padding_factor;
					zp[tid] = (s_eulers[6] * x + s_eulers[7] * y + s_eulers[8] * z) * padding_factor;
				}
				else
				{
					xp[tid] = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
					yp[tid] = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
					zp[tid] = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;
				}

				// Only consider pixels that are projected inside the sphere in output coordinates.
				//     --JZ, Nov. 26th 2018
				if ( ( xp[tid] * xp[tid] + yp[tid] * yp[tid] + zp[tid] * zp[tid] ) > max_r2_vol)
				{
					ok_for_next = 0;
				}

				if(ok_for_next)
				{
					//WAVG
					minvsigma2 = g_Minvsigma2s[pixel];
					ctf = g_ctfs[pixel];
					img_real = g_img_real[pixel];
					img_imag = g_img_imag[pixel];
					Fweight[tid] = (XFLOAT) 0.0;
					real[tid] = (XFLOAT) 0.0;
					imag[tid] = (XFLOAT) 0.0;
					XFLOAT inv_minsigma_ctf;
					if(CTF_PREMULTIPLIED)
						inv_minsigma_ctf = weight_norm_inverse * minvsigma2;
					else
						inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;

					XFLOAT temp_real, temp_imag;

					for (unsigned long itrans = 0; itrans < translation_num; itrans++)
					{
						weight = g_weights[img * translation_num + itrans];

						if (weight >= significant_weight)
						{
							weight = weight * inv_minsigma_ctf;
                                                        Fweight[tid] += weight * ctf;

							if(DATA3D)
								CpuKernels::translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
							else
								CpuKernels::translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

							real[tid] += temp_real * weight;
							imag[tid] += temp_imag * weight;
						}
					}

					//BP
					if (Fweight[tid] > (XFLOAT) 0.0)
					{

						// Only asymmetric half is stored

						if (xp[tid] < (XFLOAT) 0.0)
						{
							// Get complex conjugated hermitian symmetry pair
							xp[tid] = -xp[tid];
							yp[tid] = -yp[tid];
							zp[tid] = -zp[tid];
							imag[tid] = -imag[tid];
						}
					} // Fweight[tid] > (RFLOAT) 0.0
				}// end for if_ok_for_next
			}  // for tid

			for(int tid=0; tid<block_size; tid++)
			{
				if (Fweight[tid] > (XFLOAT) 0.0)
				{
					int x0 = floorf(xp[tid]);
					XFLOAT fx = xp[tid] - x0;
					int x1 = x0 + 1;

					int y0 = floorf(yp[tid]);
					XFLOAT fy = yp[tid] - y0;
					y0 -= mdl_inity;
					int y1 = y0 + 1;

					int z0 = floorf(zp[tid]);
					XFLOAT fz = zp[tid] - z0;
					z0 -= mdl_initz;
					int z1 = z0 + 1;

					XFLOAT mfx = (XFLOAT)1.0 - fx;
					XFLOAT mfy = (XFLOAT)1.0 - fy;
					XFLOAT mfz = (XFLOAT)1.0 - fz;

					// Locking necessary since all threads share the same back projector
					XFLOAT dd000 = mfz * mfy * mfx;
					XFLOAT dd001 = mfz * mfy *  fx;

					size_t z0MdlxMdly = (size_t)z0 * (size_t)mdl_x * (size_t)mdl_y;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y0]);

						g_model_real  [z0MdlxMdly + y0 * mdl_x + x0]+=dd000 * real[tid];
						g_model_imag  [z0MdlxMdly + y0 * mdl_x + x0]+=dd000 * imag[tid];
						g_model_weight[z0MdlxMdly + y0 * mdl_x + x0]+=dd000 * Fweight[tid];

						g_model_real  [z0MdlxMdly + y0 * mdl_x + x1]+=dd001 * real[tid];
						g_model_imag  [z0MdlxMdly + y0 * mdl_x + x1]+=dd001 * imag[tid];
						g_model_weight[z0MdlxMdly + y0 * mdl_x + x1]+=dd001 * Fweight[tid];
					}

					XFLOAT dd010 = mfz *  fy * mfx;
					XFLOAT dd011 = mfz *  fy *  fx;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y1]);

						g_model_real  [z0MdlxMdly + y1 * mdl_x + x0]+=dd010 * real[tid];
						g_model_imag  [z0MdlxMdly + y1 * mdl_x + x0]+=dd010 * imag[tid];
						g_model_weight[z0MdlxMdly + y1 * mdl_x + x0]+=dd010 * Fweight[tid];

						g_model_real  [z0MdlxMdly + y1 * mdl_x + x1]+=dd011 * real[tid];
						g_model_imag  [z0MdlxMdly + y1 * mdl_x + x1]+=dd011 * imag[tid];
						g_model_weight[z0MdlxMdly + y1 * mdl_x + x1]+=dd011 * Fweight[tid];
					}

					XFLOAT dd100 =  fz * mfy * mfx;
					XFLOAT dd101 =  fz * mfy *  fx;

					size_t z1MdlxMdly = (size_t)z1 * (size_t)mdl_x * (size_t)mdl_y;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y0]);

						g_model_real  [z1MdlxMdly + y0 * mdl_x + x0]+=dd100 * real[tid];
						g_model_imag  [z1MdlxMdly + y0 * mdl_x + x0]+=dd100 * imag[tid];
						g_model_weight[z1MdlxMdly + y0 * mdl_x + x0]+=dd100 * Fweight[tid];

						g_model_real  [z1MdlxMdly + y0 * mdl_x + x1]+=dd101 * real[tid];
						g_model_imag  [z1MdlxMdly + y0 * mdl_x + x1]+=dd101 * imag[tid];
						g_model_weight[z1MdlxMdly + y0 * mdl_x + x1]+=dd101 * Fweight[tid];
					}

					XFLOAT dd110 =  fz *  fy * mfx;
					XFLOAT dd111 =  fz *  fy *  fx;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y1]);

						g_model_real  [z1MdlxMdly + y1 * mdl_x + x0]+=dd110 * real[tid];
						g_model_imag  [z1MdlxMdly + y1 * mdl_x + x0]+=dd110 * imag[tid];
						g_model_weight[z1MdlxMdly + y1 * mdl_x + x0]+=dd110 * Fweight[tid];

						g_model_real  [z1MdlxMdly + y1 * mdl_x + x1]+=dd111 * real[tid];
						g_model_imag  [z1MdlxMdly + y1 * mdl_x + x1]+=dd111 * imag[tid];
						g_model_weight[z1MdlxMdly + y1 * mdl_x + x1]+=dd111 * Fweight[tid];
					}
				}  // Fweight[tid] > (RFLOAT) 0.0
			}  // for tid
		} // for pass
	} // for img
}

// sincos lookup table optimization. Function translatePixel calls
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions.
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty).
template < bool CTF_PREMULTIPLIED >
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backprojectRef3D(
		unsigned long imageCount,
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT* g_weights,
		XFLOAT* g_Minvsigma2s,
		XFLOAT* g_ctfs,
		unsigned long trans_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real,
		XFLOAT *g_model_imag,
		XFLOAT *g_model_weight,
		int     max_r,
		int     max_r2,
		XFLOAT   padding_factor,
		unsigned img_x,
		unsigned img_y,
		unsigned img_z,
		size_t   img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int      mdl_inity,
		int      mdl_initz,
		tbb::spin_mutex *mutexes)
{
	int img_y_half = img_y / 2;
	int img_y_half_2 = img_y_half * img_y_half;
	int img_z_half = img_z / 2;

	int max_r2_vol = max_r2 * padding_factor * padding_factor;

	// Set up the sin and cos lookup tables
	XFLOAT sin_x[trans_num][img_x], cos_x[trans_num][img_x];
	XFLOAT sin_y[trans_num][img_y], cos_y[trans_num][img_y];

	CpuKernels::computeSincosLookupTable2D(trans_num, g_trans_x, g_trans_y, img_x, img_y,
							   &sin_x[0][0], &cos_x[0][0],
							   &sin_y[0][0], &cos_y[0][0]);
	
	// Set up some other variables
	XFLOAT s_eulers[9];
	
	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;

	XFLOAT xp[img_x], yp[img_x], zp[img_x];
	XFLOAT real[img_x], imag[img_x], Fweight[img_x];

		
	for (unsigned long img=0; img<imageCount; img++) {

		for(int i = 0; i < 9; i++)
			s_eulers[i] = g_eulers[img*9+i];

		size_t mdl_x_mdl_y = (size_t)mdl_x * (size_t)mdl_y;
		size_t pixel = 0;
		for(int iy = 0; iy < img_y; iy++) {
			int y = iy;
			if (iy > img_y_half) {
				y = iy - img_y;
			}

			int y2 = y * y;
			int xmax = sqrt((XFLOAT)(img_y_half_2 - y2));  // minimize locking if possible

			memset(Fweight,0,sizeof(XFLOAT)*img_x);
			memset(real,   0,sizeof(XFLOAT)*img_x);
			memset(imag,   0,sizeof(XFLOAT)*img_x);

			for (unsigned long itrans = 0; itrans < trans_num; itrans++)
			{
				XFLOAT weight = g_weights[img * trans_num + itrans];
				if (weight < significant_weight)
					continue;

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

				#pragma omp simd
				for(int x=0; x<xmax; x++) {
					XFLOAT minvsigma2 = g_Minvsigma2s[pixel + x];
					XFLOAT ctf        = g_ctfs       [pixel + x];
					XFLOAT inv_minsigma_ctf;
					XFLOAT my_weight;
					if(CTF_PREMULTIPLIED)
					{
						inv_minsigma_ctf = weight_norm_inverse * minvsigma2;
						my_weight = weight * inv_minsigma_ctf;
						Fweight[x] += my_weight  * ctf;
					}
					else
					{
						inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;
						my_weight = weight * inv_minsigma_ctf;
						Fweight[x] += my_weight  * ctf;
					}

					XFLOAT img_real   = g_img_real   [pixel + x];
					XFLOAT img_imag   = g_img_imag   [pixel + x];

					XFLOAT ss = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
					XFLOAT cc = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

					XFLOAT shifted_real = cc * img_real - ss * img_imag;
					XFLOAT shifted_imag = cc * img_imag + ss * img_real;
					/*
					XFLOAT shifted_real, shifted_imag;
					CpuKernels::translatePixel(x, y,  g_trans_x[itrans], g_trans_y[itrans],
								   img_real, img_imag, shifted_real, shifted_imag);
					*/
					real[x] += shifted_real * my_weight;
					imag[x] += shifted_imag * my_weight;
				}
			}

			#pragma omp simd
			for(int x=0; x<xmax; x++) {
				// Get logical coordinates in the 3D map
				xp[x] = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
				yp[x] = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
				zp[x] = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;

				// Use Fweight to discard pixels that project outside the sphere
				// (pixels with Fweight <= 0 will be skipped further down)
				//     --JZ, Nov. 26th 2018
				if (xp[x]*xp[x] + yp[x]*yp[x] + zp[x]*zp[x] > max_r2_vol)
				{
					Fweight[x] = (XFLOAT) 0.0;
				}

				// Only asymmetric half is stored
				if (xp[x] < (XFLOAT) 0.0) {
					// Get complex conjugated hermitian symmetry pair
					xp[x]   = -xp[x];
					yp[x]   = -yp[x];
					zp[x]   = -zp[x];
					imag[x] = -imag[x];
				}
			}  // for x direction

			for(int x=0; x<xmax; x++){
				if (Fweight[x] <= (XFLOAT) 0.0)
					continue;

				int x0 = floorf(xp[x]);
				XFLOAT fx = xp[x] - x0;

				int y0 = floorf(yp[x]);
				XFLOAT fy = yp[x] - y0;
				y0 -= mdl_inity;

				int z0 = floorf(zp[x]);
				XFLOAT fz = zp[x] - z0;
				z0 -= mdl_initz;

				XFLOAT mfx = (XFLOAT)1.0 - fx;
				XFLOAT mfy = (XFLOAT)1.0 - fy;
				XFLOAT mfz = (XFLOAT)1.0 - fz;

				XFLOAT mfz_mfy = mfz * mfy;
				size_t z0_mdl_x_mdl_y = (size_t)z0 * mdl_x_mdl_y;
				size_t y0_mdl_x = (size_t)y0 * (size_t)mdl_x;
				size_t idx_tmp;

				XFLOAT dd000 = mfz_mfy * mfx; // mfz *  mfy *  mfx
				XFLOAT dd001 = mfz_mfy - dd000; // mfz *  mfy *  fx

				{
					tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y0]);

					idx_tmp = z0_mdl_x_mdl_y + y0_mdl_x + (size_t)x0; // z0 * mdl_x * mdl_y + y0 * mdl_x + x0;
					g_model_real  [idx_tmp]+=dd000 * real[x];
					g_model_imag  [idx_tmp]+=dd000 * imag[x];
					g_model_weight[idx_tmp]+=dd000 * Fweight[x];

					idx_tmp = idx_tmp + 1; // z0 * mdl_x * mdl_y + y0 * mdl_x + x1;
					g_model_real  [idx_tmp]+=dd001 * real[x];
					g_model_imag  [idx_tmp]+=dd001 * imag[x];
					g_model_weight[idx_tmp]+=dd001 * Fweight[x];
				}

				XFLOAT dd010 = (mfz - mfz_mfy) * mfx; // mfz *  fy *  mfx
				XFLOAT dd011 = (mfz - mfz_mfy) - dd010; // mfz *  fy *  fx

				{
					tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y0 + 1]);

					idx_tmp = z0_mdl_x_mdl_y + y0_mdl_x + (size_t)mdl_x + (size_t)x0; // z0 * mdl_x * mdl_y + y1 * mdl_x + x0;
					g_model_real  [idx_tmp]+=dd010 * real[x];
					g_model_imag  [idx_tmp]+=dd010 * imag[x];
					g_model_weight[idx_tmp]+=dd010 * Fweight[x];

					idx_tmp = idx_tmp + 1; // z0 * mdl_x * mdl_y + y1 * mdl_x + x1;
					g_model_real  [idx_tmp]+=dd011 * real[x];
					g_model_imag  [idx_tmp]+=dd011 * imag[x];
					g_model_weight[idx_tmp]+=dd011 * Fweight[x];
				}

				XFLOAT dd100 = (mfy - mfz_mfy) * mfx; // fz *  mfy *  mfx
				XFLOAT dd101 = (mfy - mfz_mfy) - dd100; // fz *  mfy *  fx
				int z1 = z0 + 1;

				{
					tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y0]);

					idx_tmp = z0_mdl_x_mdl_y + mdl_x_mdl_y + y0_mdl_x + (size_t)x0; // z1 * mdl_x * mdl_y + y0 * mdl_x + x0;
					g_model_real  [idx_tmp]+=dd100 * real[x];
					g_model_imag  [idx_tmp]+=dd100 * imag[x];
					g_model_weight[idx_tmp]+=dd100 * Fweight[x];

					idx_tmp = idx_tmp + 1; // z1 * mdl_x * mdl_y + y0 * mdl_x + x1;
					g_model_real  [idx_tmp]+=dd101 * real[x];
					g_model_imag  [idx_tmp]+=dd101 * imag[x];
					g_model_weight[idx_tmp]+=dd101 * Fweight[x];

				}

				XFLOAT dd110 = (1 - mfz - mfy + mfz_mfy) * mfx; // fz *  fy *  mfx
				XFLOAT dd111 = (1 - mfz - mfy + mfz_mfy) - dd110; // fz *  fy *  fx

				{
					tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y0 + 1]);
					idx_tmp = z0_mdl_x_mdl_y + mdl_x_mdl_y + y0_mdl_x + (size_t)mdl_x + (size_t)x0; // z1 * mdl_x * mdl_y + y1 * mdl_x + x0;
					g_model_real  [idx_tmp]+=dd110 * real[x];
					g_model_imag  [idx_tmp]+=dd110 * imag[x];
					g_model_weight[idx_tmp]+=dd110 * Fweight[x];

					idx_tmp = idx_tmp + 1; // z1 * mdl_x * mdl_y + y1 * mdl_x + x1;
					g_model_real  [idx_tmp]+=dd111 * real[x];
					g_model_imag  [idx_tmp]+=dd111 * imag[x];
					g_model_weight[idx_tmp]+=dd111 * Fweight[x];
				}
			}  // for x direction

			pixel += (size_t)img_x;
		} // for y direction
	} // for img
}

template < bool DATA3D, bool CTF_PREMULTIPLIED >
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject3D_SGD(
		unsigned long imageCount,
		int     block_size,
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
		size_t   img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz,
		tbb::spin_mutex *mutexes)
{
	int img_y_half = img_y / 2;
	int img_z_half = img_z / 2;

	int max_r2_vol = max_r2 * padding_factor * padding_factor;
	
	// Set up some variables
	XFLOAT s_eulers[9];

	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;

	//   TODO - does this really help with the call to the projector in here?
	//
	//   We collect block_size number of values before storing the results to
	//   help vectorization and control memory accesses
	XFLOAT real[block_size], imag[block_size], Fweight[block_size];
	XFLOAT ref_real[block_size], ref_imag[block_size];
	XFLOAT xp[block_size], yp[block_size], zp[block_size];
		
	for (unsigned long img=0; img<imageCount; img++) {

		for (int i = 0; i < 9; i++)
			s_eulers[i] = g_eulers[img*9+i];


		int pixel_pass_num(0);
		if(DATA3D)
			pixel_pass_num = (ceilf((float)img_xyz/(float)block_size));
		else
			pixel_pass_num = (ceilf((float)img_xyz/(float)block_size));

		for (unsigned pass = 0; pass < pixel_pass_num; pass++)   {
			memset(Fweight,0,sizeof(XFLOAT)*block_size);
//			#pragma omp simd
			for(int tid=0; tid<block_size; tid++) {
				int ok_for_next(1);  // This flag avoids continues, helping the vectorizer

				size_t pixel(0);
				if(DATA3D)
					pixel = ((size_t)pass * (size_t)block_size) + (size_t)tid;
				else
					pixel = ((size_t)pass * (size_t)block_size) + (size_t)tid;

				if (pixel >= img_xyz)
					continue;  // just doesn't make sense to proceed in this case

				int x,y,z,xy;
				XFLOAT minvsigma2, ctf, img_real, img_imag, weight;

				if(DATA3D)
				{
					z =  CpuKernels::floorfracf(pixel, (size_t)((size_t)img_x*(size_t)img_y));
					xy = pixel % (img_x*img_y);
					x =             xy  % img_x;
					y = CpuKernels::floorfracf( xy,   (size_t)img_x);

					if (z > img_z_half)
					{
						z = z - img_z;

						if(x==0)
							ok_for_next=0;
					}
				}
				else
				{
					x =             pixel % img_x;
					y = CpuKernels::floorfracf( pixel , (size_t)img_x);
				}
				if (y > img_y_half)
				{
					y = y - img_y;
				}
				// Get logical coordinates in the 3D map

				if(DATA3D)
				{
					xp[tid] = (s_eulers[0] * x + s_eulers[1] * y + s_eulers[2] * z) * padding_factor;
					yp[tid] = (s_eulers[3] * x + s_eulers[4] * y + s_eulers[5] * z) * padding_factor;
					zp[tid] = (s_eulers[6] * x + s_eulers[7] * y + s_eulers[8] * z) * padding_factor;
				}
				else
				{
					xp[tid] = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
					yp[tid] = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
					zp[tid] = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;
				}

				if (xp[tid]*xp[tid] + yp[tid]*yp[tid] + zp[tid]*zp[tid] > max_r2_vol)
				{
					ok_for_next = 0;
				}

				if(ok_for_next)
				{
					ref_real[tid] = (XFLOAT) 0.0;
					ref_imag[tid] = (XFLOAT) 0.0;

					if(DATA3D)
						projector.project3Dmodel(
							x,y,z,
							s_eulers[0], s_eulers[1], s_eulers[2],
							s_eulers[3], s_eulers[4], s_eulers[5],
							s_eulers[6], s_eulers[7], s_eulers[8],
							ref_real[tid], ref_imag[tid]);
					else
						projector.project3Dmodel(
							x,y,
							s_eulers[0], s_eulers[1],
							s_eulers[3], s_eulers[4],
							s_eulers[6], s_eulers[7],
							ref_real[tid], ref_imag[tid]);

					//WAVG
					minvsigma2 = g_Minvsigma2s[pixel];
					ctf = g_ctfs[pixel];
					img_real = g_img_real[pixel];
					img_imag = g_img_imag[pixel];
					Fweight[tid] = (XFLOAT) 0.0;
					real[tid] = (XFLOAT) 0.0;
					imag[tid] = (XFLOAT) 0.0;
					ref_real[tid] *= ctf;
					ref_imag[tid] *= ctf;
					XFLOAT inv_minsigma_ctf;

					if(CTF_PREMULTIPLIED)
						inv_minsigma_ctf = weight_norm_inverse * minvsigma2;
					else
						inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;

					XFLOAT temp_real, temp_imag;

					for (unsigned long itrans = 0; itrans < translation_num; itrans++)
					{
						weight = g_weights[img * translation_num + itrans];

						if (weight >= significant_weight)
						{
							weight = weight * inv_minsigma_ctf;
                                                        Fweight[tid] += weight * ctf;

							if(DATA3D)
								CpuKernels::translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
							else
								CpuKernels::translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

							real[tid] += (temp_real-ref_real[tid]) * weight;
							imag[tid] += (temp_imag-ref_imag[tid]) * weight;
						}
					}

					//BP
					if (Fweight[tid] > (XFLOAT) 0.0)
					{
						// Only asymmetric half is stored
						if (xp[tid] < (XFLOAT) 0.0)
						{
							// Get complex conjugated hermitian symmetry pair
							xp[tid] = -xp[tid];
							yp[tid] = -yp[tid];
							zp[tid] = -zp[tid];
							imag[tid] = -imag[tid];
						}
					} // if (Fweight[tid] > (XFLOAT) 0.0)
				} //if(ok_for_next)
			} // for tid

			for(int tid=0; tid<block_size; tid++)
			{
				if (Fweight[tid] > (XFLOAT) 0.0)
				{
					int x0 = floorf(xp[tid]);
					XFLOAT fx = xp[tid] - x0;
					int x1 = x0 + 1;

					int y0 = floorf(yp[tid]);
					XFLOAT fy = yp[tid] - y0;
					y0 -= mdl_inity;
					int y1 = y0 + 1;

					int z0 = floorf(zp[tid]);
					XFLOAT fz = zp[tid] - z0;
					z0 -= mdl_initz;
					int z1 = z0 + 1;

					XFLOAT mfx = (XFLOAT)1.0 - fx;
					XFLOAT mfy = (XFLOAT)1.0 - fy;
					XFLOAT mfz = (XFLOAT)1.0 - fz;

					XFLOAT dd000 = mfz * mfy * mfx;
					XFLOAT dd001 = mfz * mfy *  fx;

					size_t z0MdlxMdly = (size_t)z0 * (size_t)mdl_x * (size_t)mdl_y;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y0]);

						g_model_real  [z0MdlxMdly + y0 * mdl_x + x0] += dd000 * real[tid];
						g_model_imag  [z0MdlxMdly + y0 * mdl_x + x0] += dd000 * imag[tid];
						g_model_weight[z0MdlxMdly + y0 * mdl_x + x0] += dd000 * Fweight[tid];

						g_model_real  [z0MdlxMdly + y0 * mdl_x + x1] += dd001 * real[tid];
						g_model_imag  [z0MdlxMdly + y0 * mdl_x + x1] += dd001 * imag[tid];
						g_model_weight[z0MdlxMdly + y0 * mdl_x + x1] += dd001 * Fweight[tid];
					}

					XFLOAT dd010 = mfz *  fy * mfx;
					XFLOAT dd011 = mfz *  fy *  fx;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z0 * mdl_y + y1]);

						g_model_real  [z0MdlxMdly + y1 * mdl_x + x0] += dd010 * real[tid];
						g_model_imag  [z0MdlxMdly + y1 * mdl_x + x0] += dd010 * imag[tid];
						g_model_weight[z0MdlxMdly + y1 * mdl_x + x0] += dd010 * Fweight[tid];

						g_model_real  [z0MdlxMdly + y1 * mdl_x + x1] += dd011 * real[tid];
						g_model_imag  [z0MdlxMdly + y1 * mdl_x + x1] += dd011 * imag[tid];
						g_model_weight[z0MdlxMdly + y1 * mdl_x + x1] += dd011 * Fweight[tid];
					}

					XFLOAT dd100 =  fz * mfy * mfx;
					XFLOAT dd101 =  fz * mfy *  fx;

					size_t z1MdlxMdly = (size_t)z1 * (size_t)mdl_x * (size_t)mdl_y;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y0]);

						g_model_real  [z1MdlxMdly + y0 * mdl_x + x0] += dd100 * real[tid];
						g_model_imag  [z1MdlxMdly + y0 * mdl_x + x0] += dd100 * imag[tid];
						g_model_weight[z1MdlxMdly + y0 * mdl_x + x0] += dd100 * Fweight[tid];

						g_model_real  [z1MdlxMdly + y0 * mdl_x + x1] += dd101 * real[tid];
						g_model_imag  [z1MdlxMdly + y0 * mdl_x + x1] += dd101 * imag[tid];
						g_model_weight[z1MdlxMdly + y0 * mdl_x + x1] += dd101 * Fweight[tid];

					}

					XFLOAT dd110 =  fz *  fy * mfx;
					XFLOAT dd111 =  fz *  fy *  fx;

					{
						tbb::spin_mutex::scoped_lock lock(mutexes[z1 * mdl_y + y1]);

						g_model_real  [z1MdlxMdly + y1 * mdl_x + x0] += dd110 * real[tid];
						g_model_imag  [z1MdlxMdly + y1 * mdl_x + x0] += dd110 * imag[tid];
						g_model_weight[z1MdlxMdly + y1 * mdl_x + x0] += dd110 * Fweight[tid];

						g_model_real  [z1MdlxMdly + y1 * mdl_x + x1] += dd111 * real[tid];
						g_model_imag  [z1MdlxMdly + y1 * mdl_x + x1] += dd111 * imag[tid];
						g_model_weight[z1MdlxMdly + y1 * mdl_x + x1] += dd111 * Fweight[tid];
					}

				} // Fweight[tid] > (RFLOAT) 0.0
			} // for tid
		} // for pass
	} // for img
}

template < bool CTF_PREMULTIPLIED >
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject2D_SGD(
		unsigned long imageCount,
		int     block_size,
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
		int mdl_inity,
		tbb::spin_mutex *mutexes)
{
	int img_y_half = img_y / 2;

	int max_r2_out = max_r2 * padding_factor * padding_factor;

	// pre-compute sin and cos for x and y direction
	XFLOAT sin_x[translation_num][img_x], cos_x[translation_num][img_x];
	XFLOAT sin_y[translation_num][img_y], cos_y[translation_num][img_y];

	computeSincosLookupTable2D(translation_num, g_trans_x, g_trans_y,
	                           img_x, img_y,
	                           &sin_x[0][0], &cos_x[0][0],
	                           &sin_y[0][0], &cos_y[0][0]);

	// Set up some other variables
	XFLOAT s_eulers[4];

	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;

	XFLOAT xp[img_x], yp[img_x];
	XFLOAT real[img_x], imag[img_x], Fweight[img_x];

	for (unsigned long img=0; img<imageCount; img++) {

		// Copy the rotation matrix to local variables
		s_eulers[0] = g_eulers[img*9+0];
		s_eulers[1] = g_eulers[img*9+1];
		s_eulers[2] = g_eulers[img*9+3];
		s_eulers[3] = g_eulers[img*9+4];

		size_t pixel = 0;

		for(int iy = 0; iy < img_y; iy++) {
			int y = iy;
			if (iy > img_y_half) {
				y = iy - img_y;
			}

			int xmax = img_x;

			memset(Fweight,0,sizeof(XFLOAT)*img_x);
			memset(real,   0,sizeof(XFLOAT)*img_x);
			memset(imag,   0,sizeof(XFLOAT)*img_x);

			for(int x=0; x<xmax; x++) {
				//WAVG
				XFLOAT minvsigma2 = g_Minvsigma2s[pixel + x];
				XFLOAT ctf = g_ctfs[pixel + x];
				XFLOAT my_weight;

				XFLOAT ref_real = (XFLOAT) 0.0;
				XFLOAT ref_imag = (XFLOAT) 0.0;

				projector.project2Dmodel(
						x,y,
						s_eulers[0], s_eulers[1],
						s_eulers[2], s_eulers[3],
						ref_real, ref_imag);
				ref_real *= ctf;
				ref_imag *= ctf;

				for (unsigned long itrans = 0; itrans < translation_num; itrans++)
				{
					XFLOAT weight = g_weights[img * translation_num + itrans];

					if (weight < significant_weight)
						continue;

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

					if(CTF_PREMULTIPLIED) {
						my_weight = weight * weight_norm_inverse * minvsigma2;
						Fweight[x] += my_weight  * ctf;
					}
					else {
						my_weight = weight * weight_norm_inverse * ctf * minvsigma2;
						Fweight[x] += my_weight  * ctf;
					}
					/*
					CpuKernels::translatePixel(x, y,
					 * g_trans_x[itrans], g_trans_y[itrans],
					 * img_real, img_imag, temp_real, temp_imag);
					 */
					XFLOAT img_real = g_img_real[pixel + x];
					XFLOAT img_imag = g_img_imag[pixel + x];

					XFLOAT ss = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
					XFLOAT cc = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

					XFLOAT temp_real = cc * img_real - ss * img_imag;
					XFLOAT temp_imag = cc * img_imag + ss * img_real;

					real[x] += (temp_real-ref_real) * my_weight;
					imag[x] += (temp_imag-ref_imag) * my_weight;
				}  // for x
			}  // for itrans

			for(int x=0; x<xmax; x++) {
				if (Fweight[x] <= (XFLOAT) 0.0)
					continue;

				// Get logical coordinates in the 3D map
				xp[x] = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
				yp[x] = (s_eulers[2] * x + s_eulers[3] * y ) * padding_factor;

				// Only consider pixels that are projected inside the allowed circle in output coordinates.
				//     --JZ, Nov. 26th 2018
				if ( ( xp[x] * xp[x] + yp[x] * yp[x] ) > max_r2_out)
				{
					Fweight[x]= (XFLOAT) 0.0;
					continue;
				}

				// Only asymmetric half is stored
				if (xp[x] < (XFLOAT) 0.0)
				{
					// Get complex conjugated hermitian symmetry pair
					xp[x] = -xp[x];
					yp[x] = -yp[x];
					imag[x] = -imag[x];
				}
			}  // for x

			for(int x=0; x<xmax; x++) {
				if (Fweight[x] <= (XFLOAT) 0.0)
					continue;

				int x0 = floorf(xp[x]);
				XFLOAT fx = xp[x] - x0;

				int y0 = floorf(yp[x]);
				XFLOAT fy = yp[x] - y0;
				y0 -= mdl_inity;
				int y1 = y0 + 1;

				XFLOAT mfx = (XFLOAT) 1.0 - fx;
				XFLOAT mfy = (XFLOAT) 1.0 - fy;

				XFLOAT dd00 = mfy * mfx;
				XFLOAT dd01 = mfy *  fx;
				XFLOAT dd10 =  fy * mfx;
				XFLOAT dd11 =  fy *  fx;

				size_t idx_tmp;

				// Locking necessary since all threads share the same back projector
				{
					idx_tmp = (size_t)y0 * (size_t)mdl_x + (size_t)x0;

					tbb::spin_mutex::scoped_lock lock(mutexes[y0]);
					g_model_real  [idx_tmp] += dd00 * real[x];
					g_model_imag  [idx_tmp] += dd00 * imag[x];
					g_model_weight[idx_tmp] += dd00 * Fweight[x];

					idx_tmp = idx_tmp + 1;  // x1 = x0 + 1

					g_model_real  [idx_tmp] += dd01 * real[x];
					g_model_imag  [idx_tmp] += dd01 * imag[x];
					g_model_weight[idx_tmp] += dd01 * Fweight[x];
				}  // scoping for first lock

				{
					idx_tmp = (size_t)y1 * (size_t)mdl_x + (size_t)x0;

					tbb::spin_mutex::scoped_lock lock(mutexes[y1]);
					g_model_real  [idx_tmp] += dd10 * real[x];
					g_model_imag  [idx_tmp] += dd10 * imag[x];
					g_model_weight[idx_tmp] += dd10 * Fweight[x];

					idx_tmp = idx_tmp + 1;  // x1 = x0 + 1
					g_model_real  [idx_tmp] += dd11 * real[x];
					g_model_imag  [idx_tmp] += dd11 * imag[x];
					g_model_weight[idx_tmp] += dd11 * Fweight[x];
				}  // scoping for second lock
			}  // for x

			pixel += (size_t)img_x;
		} // for y
	} // for img
}


} // namespace
