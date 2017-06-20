#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <cassert>

#include "src/acc/cpu/cpu_backprojector.h"
#include "src/acc/cpu/cpu_kernels/helper.h"


namespace CpuKernels
{

#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_PREFETCH_COUNT 3
#define BP_2D_BLOCK_SIZE 128
#define BP_REF3D_BLOCK_SIZE 128
#define BP_DATA3D_BLOCK_SIZE 640

size_t Backprojector::setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int max_r, int paddingFactor)
{
	if (xdim != mdlX ||
		ydim != mdlY ||
		zdim != mdlZ ||
		inity != mdlInitY ||
		initz != mdlInitZ ||
		max_r != maxR ||
		paddingFactor != padding_factor)
	{
		clear();

		mdlX = xdim;
		mdlY = ydim;
		mdlZ = zdim;
		if (mdlZ < 1) mdlZ = 1;
		mdlXYZ = xdim*ydim*zdim;
		mdlInitY = inity;
		mdlInitZ = initz;
		maxR = max_r;
		maxR2 = max_r*max_r;
		padding_factor = paddingFactor;

		//Allocate space for model
		d_mdlReal   = (XFLOAT *)malloc(mdlXYZ * sizeof(XFLOAT));
		d_mdlImag   = (XFLOAT *)malloc(mdlXYZ * sizeof(XFLOAT));
		d_mdlWeight = (XFLOAT *)malloc(mdlXYZ * sizeof(XFLOAT));

		allocaton_size = mdlXYZ * sizeof(XFLOAT) * 3;
	}

	return allocaton_size;
}

void Backprojector::initMdl()
{
#ifdef CUDA_DEBUG
	if (mdlXYZ == 0)
	{
		printf("Model dimensions must be set with setMdlDim before call to setupMdl.");
		raise(SIGSEGV);
	}
	if (voxelCount != 0)
	{
		printf("DEBUG_ERROR: Duplicated call to model setup");
		raise(SIGSEGV);
	}
#endif

	//Initiate model with zeros
	memset(d_mdlReal,     0, mdlXYZ * sizeof(XFLOAT));
	memset(d_mdlImag,     0, mdlXYZ * sizeof(XFLOAT));
	memset(d_mdlWeight,   0, mdlXYZ * sizeof(XFLOAT));
}

void backproject2D(
		int     blockIdx_x,
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
		int mdl_inity)
{
	unsigned img = blockIdx_x;

	XFLOAT s_eulers[4];
	
	s_eulers[0] = g_eulers[img*9+0] * padding_factor;
	s_eulers[1] = g_eulers[img*9+1] * padding_factor;
	s_eulers[2] = g_eulers[img*9+3] * padding_factor;
	s_eulers[3] = g_eulers[img*9+4] * padding_factor;     
	
	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
	XFLOAT inv_minsigma_ctf;

	XFLOAT minvsigma2, ctf, img_real, img_imag, Fweight, real, imag, weight;

	int pixel_pass_num(ceilf((float)img_xy/(float)BP_2D_BLOCK_SIZE));

	for(int tid=0; tid<block_size; tid++)
	{
		for (unsigned pass = 0; pass < pixel_pass_num; pass++)
		{
			unsigned pixel = (pass * BP_2D_BLOCK_SIZE) + tid;

			if (pixel >= img_xy)
				continue;

			int x = pixel % img_x;
			int y = (int)floorf( (float)pixel / (float)img_x);

			// Don't search beyond square with side max_r
			if (y > max_r)
			{
				if (y >= img_y - max_r)
					y -= img_y;
				else
					continue;
			}

			if (x * x + y * y > max_r2)
				continue;

			//WAVG
			minvsigma2 = g_Minvsigma2s[pixel];
			ctf = g_ctfs[pixel];
			img_real = g_img_real[pixel];
			img_imag = g_img_imag[pixel];
			Fweight = (XFLOAT) 0.0;
			real = (XFLOAT) 0.0;
			imag = (XFLOAT) 0.0;
			inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;

			XFLOAT temp_real, temp_imag;
			for (unsigned long itrans = 0; itrans < translation_num; itrans++)
			{
				weight = g_weights[img * translation_num + itrans];

				if (weight >= significant_weight)
				{
					weight = weight * inv_minsigma_ctf;
					Fweight += weight * ctf;

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


				g_model_real  [y0 * mdl_x + x0]+=dd00 * real;
				g_model_imag  [y0 * mdl_x + x0]+=dd00 * imag;
				g_model_weight[y0 * mdl_x + x0]+=dd00 * Fweight;

				g_model_real  [y0 * mdl_x + x1]+=dd01 * real;
				g_model_imag  [y0 * mdl_x + x1]+=dd01 * imag;
				g_model_weight[y0 * mdl_x + x1]+=dd01 * Fweight;

				g_model_real  [y1 * mdl_x + x0]+=dd10 * real;
				g_model_imag  [y1 * mdl_x + x0]+=dd10 * imag;
				g_model_weight[y1 * mdl_x + x0]+=dd10 * Fweight;

				g_model_real  [y1 * mdl_x + x1]+=dd11 * real;
				g_model_imag  [y1 * mdl_x + x1]+=dd11 * imag;
				g_model_weight[y1 * mdl_x + x1]+=dd11 * Fweight;
			}  // Fweight > (RFLOAT) 0.0
		} // for pass
	}  // for tid
}

template < bool DATA3D >
void backproject3D(
		int     blockIdx_x,
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
		unsigned img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int mdl_inity,
		int mdl_initz)
{
	unsigned img = blockIdx_x;
	XFLOAT s_eulers[9];
	XFLOAT minvsigma2, ctf, img_real, img_imag, weight;

	 for (int i = 0; i < 9; i++)
		 s_eulers[i] = g_eulers[img*9+i];

	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;

	int pixel_pass_num(0);
	if(DATA3D)
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_DATA3D_BLOCK_SIZE));
	else
		pixel_pass_num = (ceilf((float)img_xyz/(float)BP_REF3D_BLOCK_SIZE));
	
	XFLOAT real[block_size], imag[block_size], Fweight[block_size];
	XFLOAT xp[block_size], yp[block_size], zp[block_size];
	
	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{
		memset(Fweight,0,sizeof(XFLOAT)*block_size);
		#pragma simd
		for(int tid=0; tid<block_size; tid++)
		{
			int ok_for_next(1);
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
				if (z > max_r)
				{
					if (z >= img_z - max_r)
						z = z - img_z;
					else
						ok_for_next=0;

					if(x==0)
						ok_for_next=0;	
				}
			}
			else
			{
				x =             pixel % img_x;
				y = floorfracf( pixel , img_x);
			}
			if (y > max_r)
			{
				if (y >= img_y - max_r)
					y = y - img_y;
				else
					ok_for_next=0;
			}

			if(DATA3D)
			{
				if ( ( x * x + y * y  + z * z ) > max_r2)
					ok_for_next=0;
			}
			else
			{
				if ( ( x * x + y * y ) > max_r2)
					ok_for_next=0;
			}
			//WAVG
			if(ok_for_next)
			{
				minvsigma2 = g_Minvsigma2s[pixel];
				ctf = g_ctfs[pixel];
				img_real = g_img_real[pixel];
				img_imag = g_img_imag[pixel];
				Fweight[tid] = (XFLOAT) 0.0;
				real[tid] = (XFLOAT) 0.0;
				imag[tid] = (XFLOAT) 0.0;
				XFLOAT inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;

				XFLOAT temp_real, temp_imag;
				for (unsigned long itrans = 0; itrans < translation_num; itrans++)
				{
					weight = g_weights[img * translation_num + itrans];

					if (weight >= significant_weight)
					{
						weight = weight * inv_minsigma_ctf;
						Fweight[tid] += weight * ctf;
						if(DATA3D)
							translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, temp_real, temp_imag);
						else
							translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, temp_real, temp_imag);

						real[tid] += temp_real * weight;
						imag[tid] += temp_imag * weight;
					}
				}

				//BP
				if (Fweight[tid] > (XFLOAT) 0.0)
				{
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

				   
				XFLOAT dd000 = mfz * mfy * mfx;

				g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd000 * real[tid];
				g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd000 * imag[tid];
				g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd000 * Fweight[tid];

				XFLOAT dd001 = mfz * mfy *  fx;

				g_model_real  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd001 * real[tid];
				g_model_imag  [z0 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd001 * imag[tid];
				g_model_weight[z0 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd001 * Fweight[tid];

				XFLOAT dd010 = mfz *  fy * mfx;

				g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd010 * real[tid];
				g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd010 * imag[tid];
				g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd010 * Fweight[tid];

				XFLOAT dd011 = mfz *  fy *  fx;

				g_model_real  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd011 * real[tid];
				g_model_imag  [z0 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd011 * imag[tid];
				g_model_weight[z0 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd011 * Fweight[tid];

				XFLOAT dd100 =  fz * mfy * mfx;

				g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd100 * real[tid];
				g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd100 * imag[tid];
				g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x0]+=dd100 * Fweight[tid];

				XFLOAT dd101 =  fz * mfy *  fx;

				g_model_real  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd101 * real[tid];
				g_model_imag  [z1 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd101 * imag[tid];
				g_model_weight[z1 * mdl_x * mdl_y + y0 * mdl_x + x1]+=dd101 * Fweight[tid];

				XFLOAT dd110 =  fz *  fy * mfx;

				g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd110 * real[tid];
				g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd110 * imag[tid];
				g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x0]+=dd110 * Fweight[tid];

				XFLOAT dd111 =  fz *  fy *  fx;

				g_model_real  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd111 * real[tid];
				g_model_imag  [z1 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd111 * imag[tid];
				g_model_weight[z1 * mdl_x * mdl_y + y1 * mdl_x + x1]+=dd111 * Fweight[tid];
			}  // Fweight[tid] > (RFLOAT) 0.0
		}  // for tid
	} // for pass
}

// sincos lookup table optimization. Function translatePixel calls 
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions. 
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and 
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to 
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty). 
void backprojectRef3D(
		int     img,
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
		unsigned img_xyz,
		unsigned mdl_x,
		unsigned mdl_y,
		int      mdl_inity,
		int      mdl_initz)
{
	XFLOAT s_eulers[9];
	for(int i = 0; i < 9; i++)
		s_eulers[i] = g_eulers[img*9+i];

	XFLOAT sin_x[trans_num][img_x], cos_x[trans_num][img_x];
	XFLOAT sin_y[trans_num][img_y], cos_y[trans_num][img_y];
	
	computeSincosLookupTable2D(trans_num, g_trans_x, g_trans_y, img_x, img_y,
							   &sin_x[0][0], &cos_x[0][0], 
							   &sin_y[0][0], &cos_y[0][0]);	
	
	XFLOAT weight_norm_inverse = (XFLOAT) 1.0 / weight_norm;
	
	XFLOAT xp[img_x], yp[img_x], zp[img_x];
	XFLOAT real[img_x], imag[img_x], Fweight[img_x];

	int mdl_x_mdl_y = mdl_x * mdl_y;
	int pixel = 0;
	for(int iy = 0; iy < img_y; iy++) {
		int y = iy;
		if (iy > max_r) {
			if (iy >= img_y - max_r)
				y = iy - img_y;
			else {
				pixel += img_x;
				continue;
			}
		}
				
		int y2 = y * y;
		int xmax = sqrt((XFLOAT)(max_r2 - y2));

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
							
			#pragma simd
			for(int x=0; x<xmax; x++) {
				XFLOAT minvsigma2 = g_Minvsigma2s[pixel + x];
				XFLOAT ctf        = g_ctfs       [pixel + x];
				XFLOAT inv_minsigma_ctf = weight_norm_inverse * ctf * minvsigma2;
				XFLOAT my_weight = weight * inv_minsigma_ctf;
				Fweight[x] += my_weight  * ctf;
					 
				XFLOAT img_real   = g_img_real   [pixel + x];
				XFLOAT img_imag   = g_img_imag   [pixel + x];
					 
				XFLOAT ss = trans_sin_x[x] * trans_cos_y + trans_cos_x[x] * trans_sin_y;
				XFLOAT cc = trans_cos_x[x] * trans_cos_y - trans_sin_x[x] * trans_sin_y;

				XFLOAT shifted_real = cc * img_real - ss * img_imag;
				XFLOAT shifted_imag = cc * img_imag + ss * img_real;
				/*
				XFLOAT shifted_real, shifted_imag;
				translatePixel(x, y,  g_trans_x[itrans], g_trans_y[itrans], 
							   img_real, img_imag, shifted_real, shifted_imag);
				*/
				real[x] += shifted_real * my_weight;
				imag[x] += shifted_imag * my_weight;
			}
		}

		#pragma simd
		for(int x=0; x<img_x; x++) {		
			// Get logical coordinates in the 3D map
			xp[x] = (s_eulers[0] * x + s_eulers[1] * y ) * padding_factor;
			yp[x] = (s_eulers[3] * x + s_eulers[4] * y ) * padding_factor;
			zp[x] = (s_eulers[6] * x + s_eulers[7] * y ) * padding_factor;
 
			// Only asymmetric half is stored            
			if (xp[x] < (XFLOAT) 0.0) {
				// Get complex conjugated hermitian symmetry pair
				xp[x]   = -xp[x];
				yp[x]   = -yp[x];
				zp[x]   = -zp[x];
				imag[x] = -imag[x];
			}
		}  // for x direction
 
		for(int x=0; x<img_x; x++){
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
			int z0_mdl_x_mdl_y = z0 * mdl_x_mdl_y;
			int y0_mdl_x = y0 * mdl_x;
			int idx_tmp;
				   
			XFLOAT dd000 = mfz_mfy * mfx; // mfz *  mfy *  mfx

			idx_tmp = z0_mdl_x_mdl_y + y0_mdl_x + x0; // z0 * mdl_x * mdl_y + y0 * mdl_x + x0;
			g_model_real  [idx_tmp]+=dd000 * real[x];
			g_model_imag  [idx_tmp]+=dd000 * imag[x];
			g_model_weight[idx_tmp]+=dd000 * Fweight[x];

			XFLOAT dd001 = mfz_mfy - dd000; // mfz *  mfy *  fx

			idx_tmp = idx_tmp + 1; // z0 * mdl_x * mdl_y + y0 * mdl_x + x1;
			g_model_real  [idx_tmp]+=dd001 * real[x];
			g_model_imag  [idx_tmp]+=dd001 * imag[x];
			g_model_weight[idx_tmp]+=dd001 * Fweight[x];

			XFLOAT dd010 = (mfz - mfz_mfy) * mfx; // mfz *  fy *  mfx

			idx_tmp = z0_mdl_x_mdl_y + y0_mdl_x + mdl_x + x0; // z0 * mdl_x * mdl_y + y1 * mdl_x + x0;
			g_model_real  [idx_tmp]+=dd010 * real[x];
			g_model_imag  [idx_tmp]+=dd010 * imag[x];
			g_model_weight[idx_tmp]+=dd010 * Fweight[x];

			XFLOAT dd011 = (mfz - mfz_mfy) - dd010; // mfz *  fy *  fx

			idx_tmp = idx_tmp + 1; // z0 * mdl_x * mdl_y + y1 * mdl_x + x1;
			g_model_real  [idx_tmp]+=dd011 * real[x];
			g_model_imag  [idx_tmp]+=dd011 * imag[x];
			g_model_weight[idx_tmp]+=dd011 * Fweight[x];

			XFLOAT dd100 = (mfy - mfz_mfy) * mfx; // fz *  mfy *  mfx

			idx_tmp = z0_mdl_x_mdl_y + mdl_x_mdl_y + y0_mdl_x + x0; // z1 * mdl_x * mdl_y + y0 * mdl_x + x0;
			g_model_real  [idx_tmp]+=dd100 * real[x];
			g_model_imag  [idx_tmp]+=dd100 * imag[x];
			g_model_weight[idx_tmp]+=dd100 * Fweight[x];

			XFLOAT dd101 = (mfy - mfz_mfy) - dd100; // fz *  mfy *  fx

			idx_tmp = idx_tmp + 1; // z1 * mdl_x * mdl_y + y0 * mdl_x + x1;
			g_model_real  [idx_tmp]+=dd101 * real[x];
			g_model_imag  [idx_tmp]+=dd101 * imag[x];
			g_model_weight[idx_tmp]+=dd101 * Fweight[x];

			XFLOAT dd110 = (1 - mfz - mfy + mfz_mfy) * mfx; // fz *  fy *  mfx

			idx_tmp = z0_mdl_x_mdl_y + mdl_x_mdl_y + y0_mdl_x + mdl_x + x0; // z1 * mdl_x * mdl_y + y1 * mdl_x + x0;
			g_model_real  [idx_tmp]+=dd110 * real[x];
			g_model_imag  [idx_tmp]+=dd110 * imag[x];
			g_model_weight[idx_tmp]+=dd110 * Fweight[x];

			XFLOAT dd111 = (1 - mfz - mfy + mfz_mfy) - dd110; // fz *  fy *  fx

			idx_tmp = idx_tmp + 1; // z1 * mdl_x * mdl_y + y1 * mdl_x + x1;
			g_model_real  [idx_tmp]+=dd111 * real[x];
			g_model_imag  [idx_tmp]+=dd111 * imag[x];
			g_model_weight[idx_tmp]+=dd111 * Fweight[x];
		}  // for x direction
		
		pixel += img_x;
	} // for y direction
}

void Backprojector::backproject(
		XFLOAT *d_img_real,
		XFLOAT *d_img_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT* d_weights,
		XFLOAT* d_Minvsigma2s,
		XFLOAT* d_ctfs,
		unsigned long translation_num,
		XFLOAT significant_weight,
		XFLOAT weight_norm,
		XFLOAT *d_eulers,
		int imgX,
		int imgY,
		int imgZ,
		unsigned long imageCount,
		bool data_is_3D)
{   
	if(mdlZ==1)
	{
		for(int i=0; i<(int)(imageCount); i++) {
			backproject2D(i, BP_2D_BLOCK_SIZE,
				&d_img_real[0],
				&d_img_imag[0],
				&trans_x[0],
				&trans_y[0],
				&d_weights[0],
				&d_Minvsigma2s[0],
				&d_ctfs[0],
				translation_num,
				significant_weight,
				weight_norm,
				d_eulers,
				d_mdlReal,
				d_mdlImag,
				d_mdlWeight,
				maxR,
				maxR2,
				padding_factor,
				imgX,
				imgY,
				imgX*imgY,
				mdlX,
				mdlInitY);
		}
	}
	else
	{
		if(data_is_3D)
			for(int i=0; i<(int)(imageCount); i++) {
				backproject3D<true>(i, BP_DATA3D_BLOCK_SIZE,
				&d_img_real[0],
				&d_img_imag[0],
				&trans_x[0],
				&trans_y[0],
				&trans_z[0],
				&d_weights[0],
				&d_Minvsigma2s[0],
				&d_ctfs[0],
				translation_num,
				significant_weight,
				weight_norm,
				d_eulers,
				d_mdlReal,
				d_mdlImag,
				d_mdlWeight,
				maxR,
				maxR2,
				padding_factor,
				imgX,
				imgY,
				imgZ,
				imgX*imgY*imgZ,
				mdlX,
				mdlY,
				mdlInitY,
				mdlInitZ);
			}
		else {
			for(int i=0; i<(int)(imageCount); i++) {
				backprojectRef3D(i,
				&d_img_real[0],
				&d_img_imag[0],
				&trans_x[0],
				&trans_y[0],
				&d_weights[0],
				&d_Minvsigma2s[0],
				&d_ctfs[0],
				translation_num,
				significant_weight,
				weight_norm,
				d_eulers,
				d_mdlReal,
				d_mdlImag,
				d_mdlWeight,
				maxR,
				maxR2,
				padding_factor,
				imgX,
				imgY,
				imgZ,
				imgX*imgY*imgZ,
				mdlX,
				mdlY,
				mdlInitY,
				mdlInitZ);
			}
		}		
	}
}


void Backprojector::getMdlData(XFLOAT *r, XFLOAT *i, XFLOAT * w)
{
	memcpy(r, d_mdlReal,   mdlXYZ * sizeof(XFLOAT));
	memcpy(i, d_mdlImag,   mdlXYZ * sizeof(XFLOAT));
	memcpy(w, d_mdlWeight, mdlXYZ * sizeof(XFLOAT));
}

void Backprojector::clear()
{
	mdlX = 0;
	mdlY = 0;
	mdlZ = 0;
	mdlXYZ = 0;
	mdlInitY = 0;
	mdlInitZ = 0;
	maxR = 0;
	maxR2 = 0;
	padding_factor = 0;
	allocaton_size = 0;

	if (d_mdlReal != NULL)
	{
		free(d_mdlReal);
		free(d_mdlImag);
		free(d_mdlWeight);

		d_mdlReal = d_mdlImag = d_mdlWeight = NULL;
	}
}

Backprojector::~Backprojector()
{
	clear();
}

} // end of namespace CpuKernels
