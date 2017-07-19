#ifndef CUDA_HELPER_KERNELS_CUH_
#define CUDA_HELPER_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projectorkernel_impl.h"

#ifdef ACC_DOUBLE_PRECISION
#define FAILSAFE_PRIOR_MIN_LIM 1e-300
#else
#define FAILSAFE_PRIOR_MIN_LIM 1e-30
#endif

template<bool failsafe,typename weights_t>
__global__ void cuda_kernel_exponentiate_weights_coarse(
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		weights_t *g_Mweight,
		XFLOAT avg_diff2,
		XFLOAT min_diff2,
		int nr_coarse_orient,
		int nr_coarse_trans)
{
	// blockid
	int bid  = blockIdx.x;
	int cid  = blockIdx.y;
	//threadid
	int tid = threadIdx.x;

	int pos, iorient = bid*SUMW_BLOCK_SIZE+tid;

	weights_t weight;
	if(iorient<nr_coarse_orient)
	{
		for (int itrans=0; itrans<nr_coarse_trans; itrans++)
		{
			pos = cid * nr_coarse_orient * nr_coarse_trans + iorient * nr_coarse_trans + itrans;
			weights_t diff2 = g_Mweight[pos];
			if( diff2 < min_diff2 ) //TODO Might be slow (divergent threads)
				weight = (weights_t)0.0;
			else
			{
				diff2 -= avg_diff2;
				weight = g_pdf_orientation[iorient] * g_pdf_offset[itrans];          	// Same for all threads - TODO: should be done once for all trans through warp-parallel execution

				if (failsafe && weight < FAILSAFE_PRIOR_MIN_LIM) //Prevent zero priors in fail-safe mode
					weight = FAILSAFE_PRIOR_MIN_LIM;

				// next line because of numerical precision of exp-function
#ifdef ACC_DOUBLE_PRECISION
				if (diff2 > 700.)
					weight = 0.;
				else
					weight *= exp(-diff2);
#else
				if (diff2 > 88.)
					weight = 0.;
				else
					weight *= expf(-diff2);
#endif
				// TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn
			}

			// Store the weight
			g_Mweight[pos] = weight; // TODO put in shared mem
		}
	}
}

template<bool DATA3D>
__global__ void cuda_kernel_collect2jobs(	XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
											XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
											XFLOAT *g_oo_otrans_z,          // otrans-size -> make const
											XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
											XFLOAT *g_i_weights,
											XFLOAT op_significant_weight,    // TODO Put in const
											XFLOAT op_sum_weight,            // TODO Put in const
											int   coarse_trans,
											int   oversamples_trans,
											int   oversamples_orient,
											int   oversamples,
											bool  do_ignore_pdf_direction,
											XFLOAT *g_o_weights,
											XFLOAT *g_thr_wsum_prior_offsetx_class,
											XFLOAT *g_thr_wsum_prior_offsety_class,
											XFLOAT *g_thr_wsum_prior_offsetz_class,
											XFLOAT *g_thr_wsum_sigma2_offset,
											unsigned long *d_rot_idx,
											unsigned long *d_trans_idx,
											unsigned long *d_job_idx,
											unsigned long *d_job_num
											)
{
	// blockid
	int bid = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	extern __shared__ XFLOAT buffer[];

	XFLOAT * s_o_weights 					= &buffer[                0];
	XFLOAT * s_thr_wsum_sigma2_offset       = &buffer[  SUMW_BLOCK_SIZE];
	XFLOAT * s_thr_wsum_prior_offsetx_class = &buffer[2*SUMW_BLOCK_SIZE];
	XFLOAT * s_thr_wsum_prior_offsety_class = &buffer[3*SUMW_BLOCK_SIZE];
	XFLOAT * s_thr_wsum_prior_offsetz_class(0);

	if(DATA3D)
		s_thr_wsum_prior_offsetz_class = &buffer[4*SUMW_BLOCK_SIZE];

	s_o_weights[tid]                    	= (XFLOAT)0.0;
	s_thr_wsum_sigma2_offset[tid]       	= (XFLOAT)0.0;

	s_thr_wsum_prior_offsetx_class[tid] 	= (XFLOAT)0.0;
	s_thr_wsum_prior_offsety_class[tid] 	= (XFLOAT)0.0;
	if(DATA3D)
		s_thr_wsum_prior_offsety_class[tid] = (XFLOAT)0.0;

	long int pos = d_job_idx[bid];
    int job_size = d_job_num[bid];
	pos += tid;	   					// pos is updated to be thread-resolved

    int pass_num = ceilfracf(job_size,SUMW_BLOCK_SIZE);
    __syncthreads();
    for (int pass = 0; pass < pass_num; pass++, pos+=SUMW_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
    {
    	if ((pass*SUMW_BLOCK_SIZE+tid)<job_size) // if there is a translation that needs to be done still for this thread
    	{
			// index of comparison
			long int iy = d_trans_idx[pos];              // ...and its own trans...

			XFLOAT weight = g_i_weights[pos];
			if( weight >= op_significant_weight ) //TODO Might be slow (divergent threads)
				weight /= op_sum_weight;
			else
				weight = (XFLOAT)0.0;

			s_o_weights[tid] 					+= weight;
			s_thr_wsum_sigma2_offset[tid]       += weight * g_myp_oo_otrans_x2y2z2[iy];
			s_thr_wsum_prior_offsetx_class[tid] += weight *          g_oo_otrans_x[iy];
			s_thr_wsum_prior_offsety_class[tid] += weight *          g_oo_otrans_y[iy];
			if(DATA3D)
				s_thr_wsum_prior_offsetz_class[tid] += weight *          g_oo_otrans_z[iy];

    	}
    }
    __syncthreads();
    // Reduction of all treanslations this orientation
	for(int j=(SUMW_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
		{
			s_o_weights[tid]                    += s_o_weights[tid+j];
			s_thr_wsum_sigma2_offset[tid]       += s_thr_wsum_sigma2_offset[tid+j];
			s_thr_wsum_prior_offsetx_class[tid] += s_thr_wsum_prior_offsetx_class[tid+j];
			s_thr_wsum_prior_offsety_class[tid] += s_thr_wsum_prior_offsety_class[tid+j];
			if(DATA3D)
				s_thr_wsum_prior_offsetz_class[tid] += s_thr_wsum_prior_offsetz_class[tid+j];
		}
		__syncthreads();
	}
	g_o_weights[bid]			        = s_o_weights[0];
	g_thr_wsum_sigma2_offset[bid]       = s_thr_wsum_sigma2_offset[0];
	g_thr_wsum_prior_offsetx_class[bid] = s_thr_wsum_prior_offsetx_class[0];
	g_thr_wsum_prior_offsety_class[bid] = s_thr_wsum_prior_offsety_class[0];
	if(DATA3D)
		g_thr_wsum_prior_offsetz_class[bid] = s_thr_wsum_prior_offsetz_class[0];
}

__global__ void cuda_kernel_exponentiate_weights_fine(    XFLOAT *g_pdf_orientation,
											  XFLOAT *g_pdf_offset,
											  XFLOAT *g_weights,
											  XFLOAT avg_diff2,
											  int oversamples_orient,
											  int oversamples_trans,
									     	  unsigned long *d_rot_id,
											  unsigned long *d_trans_idx,
											  unsigned long *d_job_idx,
											  unsigned long *d_job_num,
									     	  long int job_num);

__global__ void cuda_kernel_softMaskOutsideMap(	XFLOAT *vol,
												long int vol_size,
												long int xdim,
												long int ydim,
												long int zdim,
												long int xinit,
												long int yinit,
												long int zinit,
												bool do_Mnoise,
												XFLOAT radius,
												XFLOAT radius_p,
												XFLOAT cosine_width	);

__global__ void cuda_kernel_softMaskBackgroundValue(	XFLOAT *vol,
														long int vol_size,
														long int xdim,
														long int ydim,
														long int zdim,
														long int xinit,
														long int yinit,
														long int zinit,
														bool do_Mnoise,
														XFLOAT radius,
														XFLOAT radius_p,
														XFLOAT cosine_width,
														XFLOAT *g_sum,
														XFLOAT *g_sum_bg);

__global__ void cuda_kernel_cosineFilter(	XFLOAT *vol,
											long int vol_size,
											long int xdim,
											long int ydim,
											long int zdim,
											long int xinit,
											long int yinit,
											long int zinit,
											bool do_Mnoise,
											XFLOAT radius,
											XFLOAT radius_p,
											XFLOAT cosine_width,
											XFLOAT sum_bg_total);

//----------------------------------------------------------------------------
namespace CudaKernels
{

template <typename T>
__global__ void cuda_kernel_translate2D(	T * g_image_in,
											T * g_image_out,
											int image_size,
											int xdim,
											int ydim,
											int dx,
											int dy)
{
	int tid = threadIdx.x;
	int bid =  blockIdx.x;

	int x,y,xp,yp;
	int pixel=tid + bid*BLOCK_SIZE;
	int new_pixel;

	if(pixel<image_size)
	{
		x = pixel % xdim;
		y = (pixel-x) / (xdim);

		xp = x + dx;
		yp = y + dy;

		if( yp>=0 && xp>=0 && yp<ydim && xp<xdim)
		{
			new_pixel = yp*xdim + xp;
			if(new_pixel>=0 && new_pixel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_pixel] = g_image_in[pixel];
		}
	}
}

template <typename T>
__global__ void cuda_kernel_translate3D(	T * g_image_in,
											T * g_image_out,
											int image_size,
											int xdim,
											int ydim,
											int zdim,
											int dx,
											int dy,
											int dz)
{
	int tid = threadIdx.x;
	int bid =  blockIdx.x;

	int x,y,z,xp,yp,zp,xy;
	int voxel=tid + bid*BLOCK_SIZE;
	int new_voxel;

	int xydim = xdim*ydim;

	if(voxel<image_size)
	{
		z =  voxel / xydim;
		zp = z + dz;

		xy = voxel % xydim;
		y =  xy / xdim;
		yp = y + dy;

		x =  xy % xdim;
		xp = x + dx;

		if( zp>=0 && yp>=0 && xp>=0 && zp<zdim && yp<ydim && xp<xdim)
		{
			new_voxel = zp*xydim +  yp*xdim + xp;
			if(new_voxel>=0 && new_voxel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_voxel] = g_image_in[voxel];
		}
	}
}
} // namespace CudaKernels

//----------------------------------------------------------------------------
//__global__ void cuda_kernel_selfTranslate2D(	XFLOAT * g_image_in,
//												XFLOAT * g_image_out,
//												int image_size,
//												int xdim,
//												int ydim, //not used
//												int dx,
//												int dy);
//
//__global__ void cuda_kernel_selfTranslate3D(	XFLOAT * g_image_in,
//												XFLOAT * g_image_out,
//												int image_size,
//												int xdim,
//												int ydim,
//												int zdim, //not used
//												int dx,
//												int dy,
//												int dz);
//----------------------------------------------------------------------------
//__global__ void cuda_kernel_powerClass2D(	ACCCOMPLEX * g_image,
//											XFLOAT * g_spectrum,
//											int image_size,
//											int spectrum_size,
//											int xdim,
//											int ydim,
//											int res_limit,
//											XFLOAT * g_highres_Xi2);
//
//__global__ void cuda_kernel_powerClass3D(	ACCCOMPLEX * g_image,
//											XFLOAT * g_spectrum,
//											int image_size,
//											int spectrum_size,
//											int xdim,
//											int ydim,
//											int zdim,
//											int res_limit,
//											XFLOAT * g_highres_Xi2);

//----------------------------------------------------------------------------
__global__ void cuda_kernel_centerFFT_2D(XFLOAT *img_in,
										 int image_size,
										 int xdim,
										 int ydim,
										 int xshift,
										 int yshift);

__global__ void cuda_kernel_centerFFT_3D(XFLOAT *img_in,
										 int image_size,
										 int xdim,
										 int ydim,
										 int zdim,
										 int xshift,
										 int yshift,
										 int zshift);
//----------------------------------------------------------------------------
__global__ void cuda_kernel_probRatio(  XFLOAT *d_Mccf,
										XFLOAT *d_Mpsi,
										XFLOAT *d_Maux,
										XFLOAT *d_Mmean,
										XFLOAT *d_Mstddev,
										int image_size,
										XFLOAT normfft,
										XFLOAT sum_ref_under_circ_mask,
										XFLOAT sum_ref2_under_circ_mask,
										XFLOAT expected_Pratio,
										int NpsiThisBatch,
										int startPsi,
										int totalPsis);

__global__ void cuda_kernel_rotateOnly(   ACCCOMPLEX *d_Faux,
						  	  	  	  	  XFLOAT psi,
						  	  			  AccProjectorKernel projector,
						  	  			  int startPsi);

__global__ void cuda_kernel_rotateAndCtf( ACCCOMPLEX *d_Faux,
						  	  	  	  	  XFLOAT *d_ctf,
						  	  	  	  	  XFLOAT psi,
						  	  			  AccProjectorKernel projector,
						  	  			  int startPsi = 0);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A
 */
__global__ void cuda_kernel_convol_A(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 int image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A, writes to C
 */
__global__ void cuda_kernel_convol_A(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 ACCCOMPLEX *d_C,
									 	 int image_size);

/*
 * Multiplies many complex arrays A (in-place) by a single B, pixel-by-pixel, after conjugating A
 */
__global__ void cuda_kernel_batch_convol_A(	 ACCCOMPLEX *d_A,
									 	 	 ACCCOMPLEX *d_B,
									 	 	 int image_size);

/*
* Multiplies many complex arrays A (not in-place) by a single B, pixel-by-pixel, after conjugating A
*/
__global__ void cuda_kernel_batch_convol_A( ACCCOMPLEX *d_A,
								ACCCOMPLEX *d_B,
								ACCCOMPLEX *d_C,
								int image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B
 */
__global__ void cuda_kernel_convol_B(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 int image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B, writes to C
 */
__global__ void cuda_kernel_convol_B(	 ACCCOMPLEX *d_A,
									 	 ACCCOMPLEX *d_B,
									 	 ACCCOMPLEX *d_C,
									 	 int image_size);
/*
 * Multiplies many complex arrays A (in-place) by a single one B, pixel-by-pixel, after conjugating B
 */
__global__ void cuda_kernel_batch_convol_B(	 ACCCOMPLEX *d_A,
									 	 	 ACCCOMPLEX *d_B,
									 	 	 int image_size);

/*
 * Multiplies scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi( T *A,
								   T *OUT,
								   T S,
		  	  	  	  	  	  	   int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*S;
}

namespace CudaKernels
{
/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi(
		T *A,
		T S,
		int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*S;
}
}
/*
 * Multiplies scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi( T *A,
								   T *B,
								   T *OUT,
								   T S,
		  	  	  	  	  	  	   int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*B[pixel]*S;
}

__global__ void cuda_kernel_finalizeMstddev( XFLOAT *Mstddev,
											 XFLOAT *aux,
											 XFLOAT S,
											 int image_size);

/*
 * In place squares array in place
 *
 *  A[i] = A[i]*A[i]
 */
__global__ void cuda_kernel_square(
		XFLOAT *A,
		int image_size);

/*
 * Casts on device so we can copy_to_host directly into a multidimarray.
 */
template <typename T1, typename T2 >
__global__ void cuda_kernel_cast(
		T1 *IN,
		T2 *OUT,
		int size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<size)
		OUT[pixel] = IN[pixel];
}

template<bool do_highpass>
__global__ void cuda_kernel_frequencyPass(
										ACCCOMPLEX *A,
										long int ori_size,
										size_t Xdim,
										size_t Ydim,
										size_t Zdim,
										XFLOAT edge_low,
										XFLOAT edge_width,
										XFLOAT edge_high,
										XFLOAT angpix,
										int image_size)
{
	int texel = threadIdx.x + blockIdx.x*BLOCK_SIZE;

	int z = texel / (Xdim*Ydim);
	int xy = (texel - z*Xdim*Ydim);
	int y = xy / Xdim;

	int xp = xy - y*Xdim;

	int zp = ( z<Xdim ? z : z-Zdim );
	int yp = ( y<Xdim ? y : y-Ydim );

	int r2 = xp*xp + yp*yp + zp*zp;

	RFLOAT res;
	if(texel<image_size)
	{
		res = sqrt((RFLOAT)r2)/(RFLOAT)ori_size;

		if(do_highpass) //highpass
		{
			if (res < edge_low) //highpass => lows are dead
			{
				A[texel].x = 0.;
				A[texel].y = 0.;
			}
			else if (res < edge_high) //highpass => medium lows are almost dead
			{
				XFLOAT mul = 0.5 - 0.5 * cos( PI * (res-edge_low)/edge_width);
				A[texel].x *= mul;
				A[texel].y *= mul;
			}
		}
		else //lowpass
		{
			if (res > edge_high) //lowpass => highs are dead
			{
				A[texel].x = 0.;
				A[texel].y = 0.;
			}
			else if (res > edge_low) //lowpass => medium highs are almost dead
			{
				XFLOAT mul = 0.5 + 0.5 * cos( PI * (res-edge_low)/edge_width);
				A[texel].x *= mul;
				A[texel].y *= mul;
			}
		}
	}
}

template<bool DATA3D>
__global__ void cuda_kernel_powerClass(		ACCCOMPLEX * g_image,
											XFLOAT * g_spectrum,
											int image_size,
											int spectrum_size,
											int xdim,
											int ydim,
											int zdim,
											int res_limit,
											XFLOAT * g_highres_Xi2)
{
	int tid = threadIdx.x;
	int bid =  blockIdx.x;

	XFLOAT normFaux;
	__shared__ XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];
	s_highres_Xi2[tid] = (XFLOAT)0.;

	int x,y,xy,d;
	int xydim = xdim*ydim;
	int voxel=tid + bid*POWERCLASS_BLOCK_SIZE;
	bool coords_in_range(true);

	if(voxel<image_size)
	{
		if(DATA3D)
		{
			int z =  voxel / xydim;
			xy = voxel % xydim;
			y =  xy / xdim;
			x =  xy % xdim;

			y = ((y<xdim) ? y : y-ydim);
			z = ((z<xdim) ? z : z-zdim);

			d  = (x*x + y*y + z*z);
			coords_in_range = !(x==0 && y<0.f && z<0.f);
		}
		else
		{
			x = voxel % xdim;
			y = (voxel-x) / (xdim);

			y = ((y<xdim) ? y : y-ydim);
			d  = (x*x + y*y);
			coords_in_range = !(x==0 && y<0.f);
		}

#if defined(ACC_DOUBLE_PRECISION)
		int ires = __double2int_rn(sqrt((XFLOAT)d));
#else
		int ires = __float2int_rn(sqrtf((XFLOAT)d));
#endif
		if((ires>0.f) && (ires<spectrum_size) && coords_in_range)
		{
			normFaux = g_image[voxel].x*g_image[voxel].x + g_image[voxel].y*g_image[voxel].y;
			cuda_atomic_add(&g_spectrum[ires], normFaux);
			if(ires>=res_limit)
				s_highres_Xi2[tid] = normFaux;
		}
	}

	// Reduce the higres_Xi2-values for all threads. (I tried a straight atomic-write: for 128 threads it was ~3x slower)
	__syncthreads();
	for(int j=(POWERCLASS_BLOCK_SIZE/2); j>0.f; j/=2)
	{
		if(tid<j)
			s_highres_Xi2[tid] += s_highres_Xi2[tid+j];
		__syncthreads();
	}
	if(tid==0)
		cuda_atomic_add(&g_highres_Xi2[0], s_highres_Xi2[0]);
}

template<bool invert>
__global__ void cuda_kernel_make_eulers_2D(
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned orientation_num);

template<bool invert,bool perturb>
__global__ void cuda_kernel_make_eulers_3D(
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned orientation_num,
		XFLOAT *R);

#define INIT_VALUE_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_init_complex_value(
		T *data,
		XFLOAT value,
		size_t size)
{
	size_t idx = blockIdx.x * INIT_VALUE_BLOCK_SIZE + threadIdx.x;
	if (idx < size)
	{
		data[idx].x = value;
		data[idx].y = value;
	}
}

template< typename T>
__global__ void cuda_kernel_init_value(
		T *data,
		T value,
		size_t size)
{
	size_t idx = blockIdx.x * INIT_VALUE_BLOCK_SIZE + threadIdx.x;
	if (idx < size)
		data[idx] = value;
}

#define WEIGHT_MAP_BLOCK_SIZE 512
__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num,
        int block_size
		);

#define OVER_THRESHOLD_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_array_over_threshold(
		T *data,
		bool *passed,
		T threshold,
		size_t size)
{
	size_t idx = blockIdx.x * OVER_THRESHOLD_BLOCK_SIZE + threadIdx.x;
	if (idx < size)
	{
		if (data[idx] >= threshold)
			passed[idx] = true;
		else
			passed[idx] = false;
	}
}

#define FIND_IN_CUMULATIVE_BLOCK_SIZE 512
template< typename T>
__global__ void cuda_kernel_find_threshold_idx_in_cumulative(
		T *data,
		T threshold,
		size_t size_m1, //data size minus 1
		size_t *idx)
{
	size_t i = blockIdx.x * FIND_IN_CUMULATIVE_BLOCK_SIZE + threadIdx.x;
	if (i < size_m1 && data[i] <= threshold && threshold < data[i+1])
		idx[0] = i+1;
}

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
		XFLOAT *g_in_real,
		XFLOAT *g_in_imag,
		XFLOAT *g_out_real,
		XFLOAT *g_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, unsigned iYX, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ, unsigned oYX, //Output dimensions
		unsigned max_idx,
		unsigned max_r2 = 0
		)
{
	unsigned n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
	long int image_offset = oX*oY*oZ*blockIdx.y;
	if (n >= max_idx) return;

	int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	g_out_real[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_real[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
	g_out_imag[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp + image_offset] = g_in_imag[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp + image_offset];
}

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
		ACCCOMPLEX *g_in,
		ACCCOMPLEX *g_out,
		size_t iX, size_t iY, size_t iZ, size_t iYX, //Input dimensions
		size_t oX, size_t oY, size_t oZ, size_t oYX, //Output dimensions
		size_t max_idx,
		size_t max_r2 = 0
		)
{
	size_t n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
	size_t oOFF = oX*oY*oZ*blockIdx.y;
	size_t iOFF = iX*iY*iZ*blockIdx.y;
	if (n >= max_idx) return;

	long int k, i, kp, ip, jp;

	if (check_max_r2)
	{
		k = n / (iX * iY);
		i = (n % (iX * iY)) / iX;

		kp = k < iX ? k : k - iZ;
		ip = i < iX ? i : i - iY;
		jp = n % iX;

		if (kp*kp + ip*ip + jp*jp > max_r2)
			return;
	}
	else
	{
		k = n / (oX * oY);
		i = (n % (oX * oY)) / oX;

		kp = k < oX ? k : k - oZ;
		ip = i < oX ? i : i - oY;
		jp = n % oX;
	}

	long int  in_idx = (kp < 0 ? kp + iZ : kp) * iYX + (ip < 0 ? ip + iY : ip)*iX + jp;
	long int out_idx = (kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp;
	g_out[out_idx + oOFF] =  g_in[in_idx + iOFF];
}

#endif /* CUDA_HELPER_KERNELS_CUH_ */
