#ifndef HELPER_KERNELS_H_
#define HELPER_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include "src/macros.h"
#include "src/acc/cpu/cpu_settings.h"
#include "src/acc/cpu/cpu_kernels/cpu_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projectorkernel_impl.h"

namespace CpuKernels
{

template<typename T>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void weights_exponent_coarse(
		T *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		T *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		T *g_weights,
		T g_min_diff2,
		unsigned long  nr_coarse_orient,
		unsigned long  nr_coarse_trans,
		size_t max_idx)
{
	for (size_t idx = 0; idx < max_idx; idx++)
	{
		unsigned long  itrans = idx % nr_coarse_trans;
		unsigned long  iorient = (idx - itrans) / nr_coarse_trans;

		T diff2 = g_weights[idx];
		if( diff2 < g_min_diff2 || g_pdf_orientation_zeros[iorient] || g_pdf_offset_zeros[itrans])
	// TODO - replace with lowest() when C++11 is supported
			g_weights[idx] = -std::numeric_limits<T>::max(); //large negative number
		else
			g_weights[idx] = g_pdf_orientation[iorient] + g_pdf_offset[itrans] + g_min_diff2 - diff2;
	}
}


template<typename T>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void exponentiate(
		T *g_array,
		T add,
		size_t size)
{
	for (size_t idx = 0; idx < size; idx++)
	{
		T a = g_array[idx] + add;
#ifdef ACC_DOUBLE_PRECISION
		if (a < -700.)
			g_array[idx] = 0.;
		else
			g_array[idx] = exp(a);
#else
		if (a < -88.f)
			g_array[idx] = 0.f;
		else
			g_array[idx] = expf(a);
#endif
	}
}

template<bool DATA3D>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void collect2jobs(  int     grid_size,
					int		block_size,
					XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
					XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
					XFLOAT *g_oo_otrans_z,          // otrans-size -> make const
					XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
					XFLOAT *g_i_weights,
					XFLOAT op_significant_weight,    // TODO Put in const
					XFLOAT op_sum_weight,            // TODO Put in const
					unsigned long   coarse_trans,
					unsigned long   oversamples_trans,
					unsigned long   oversamples_orient,
					unsigned long   oversamples,
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
	// block id
	for (int bid=0; bid < grid_size; bid++) {

		XFLOAT s_o_weights[block_size];
		XFLOAT s_thr_wsum_sigma2_offset[block_size];;
		XFLOAT s_thr_wsum_prior_offsetx_class[block_size];
		XFLOAT s_thr_wsum_prior_offsety_class[block_size];
		XFLOAT s_thr_wsum_prior_offsetz_class[block_size];

		unsigned long pos = d_job_idx[bid];
		unsigned long job_size = d_job_num[bid];

		int pass_num = ceilfracf(job_size,block_size);

		for(int tid=0; tid<block_size; tid++) {
			s_o_weights[tid]                    	= (XFLOAT)0.0;
			s_thr_wsum_sigma2_offset[tid]       	= (XFLOAT)0.0;
			s_thr_wsum_prior_offsetx_class[tid] 	= (XFLOAT)0.0;
			s_thr_wsum_prior_offsety_class[tid] 	= (XFLOAT)0.0;
			if(DATA3D)
				s_thr_wsum_prior_offsety_class[tid] = (XFLOAT)0.0;
		}


		for (int pass = 0; pass < pass_num; pass++, pos+=block_size) // loop the available warps enough to complete all translations for this orientation
		{
			for(int tid=0; tid<block_size; tid++) {
				if ((pass*block_size+tid)<job_size) // if there is a translation that needs to be done still for this thread
				{
					// index of comparison
					long int iy = d_trans_idx[pos+tid];              // ...and its own trans...

					XFLOAT weight = g_i_weights[pos+tid];
					if( weight >= op_significant_weight ) //TODO Might be slow (divergent threads)
						weight /= op_sum_weight;
					else
						weight = (XFLOAT)0.0;

					s_o_weights[tid] += weight;
					s_thr_wsum_prior_offsetx_class[tid] += weight *          g_oo_otrans_x[iy];
					s_thr_wsum_prior_offsety_class[tid] += weight *          g_oo_otrans_y[iy];
					s_thr_wsum_sigma2_offset[tid]       += weight * g_myp_oo_otrans_x2y2z2[iy];
				}
			} 
		}

		for(int tid=1; tid<block_size; tid++)
		{
			s_o_weights[0]                    += s_o_weights[tid];
			s_thr_wsum_sigma2_offset[0]       += s_thr_wsum_sigma2_offset[tid];
			s_thr_wsum_prior_offsetx_class[0] += s_thr_wsum_prior_offsetx_class[tid];
			s_thr_wsum_prior_offsety_class[0] += s_thr_wsum_prior_offsety_class[tid];
			if(DATA3D)
				s_thr_wsum_prior_offsetz_class[0] += s_thr_wsum_prior_offsetz_class[tid];
		}
		g_o_weights[bid]			        = s_o_weights[0];
		g_thr_wsum_sigma2_offset[bid]       = s_thr_wsum_sigma2_offset[0];
		g_thr_wsum_prior_offsetx_class[bid] = s_thr_wsum_prior_offsetx_class[0];
		g_thr_wsum_prior_offsety_class[bid] = s_thr_wsum_prior_offsety_class[0];
		if(DATA3D)
			g_thr_wsum_prior_offsetz_class[bid] = s_thr_wsum_prior_offsetz_class[0];
	} // for bid
}

void exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		bool *g_pdf_orientation_zeros,
		XFLOAT *g_pdf_offset,
		bool *g_pdf_offset_zeros,
		XFLOAT *g_weights,
		XFLOAT min_diff2,
		unsigned long oversamples_orient,
		unsigned long oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num);

void RNDnormalDitributionComplexWithPowerModulation2D(ACCCOMPLEX* Image, size_t xdim, XFLOAT *spectra);
void RNDnormalDitributionComplexWithPowerModulation3D(ACCCOMPLEX* Image, size_t xdim, size_t ydim, XFLOAT *spectra);

void softMaskBackgroundValue(	int      block_dim,
                                int      block_size,
                                XFLOAT  *vol,
								long int vol_size,
								long int xdim,
								long int ydim,
								long int zdim,
								long int xinit,
								long int yinit,
								long int zinit,
								XFLOAT   radius,
								XFLOAT   radius_p,
								XFLOAT   cosine_width,
								XFLOAT  *g_sum,
								XFLOAT  *g_sum_bg);

void cosineFilter(	int      block_dim,
                    int      block_size,
					XFLOAT *vol,
					long int vol_size,
					long int xdim,
					long int ydim,
					long int zdim,
					long int xinit,
					long int yinit,
					long int zinit,
					bool do_Mnoise,
					XFLOAT *noise,
					XFLOAT radius,
					XFLOAT radius_p,
					XFLOAT cosine_width,
					XFLOAT sum_bg_total);
								
//----------------------------------------------------------------------------

template <typename T>
void cpu_translate2D(T * g_image_in,
					T *      g_image_out,
					size_t   image_size,
					int      xdim,
					int      ydim, //not used
					int      dx,
					int      dy);

template <typename T>
void cpu_translate3D(T * g_image_in,
					T *      g_image_out,
					size_t   image_size,
					int      xdim,
					int      ydim,
					int      zdim, //not used
					int      dx,
					int      dy,
					int      dz);

//----------------------------------------------------------------------------
template <typename T>
void centerFFT_2D(  int		batch_size,
					size_t		pixel_start,
					size_t		pixel_end,
					T		*img_in,
					size_t	image_size,
					int		xdim,
					int		ydim,
					int		xshift,
					int		yshift);

template <typename T>
void centerFFT_3D(  int		batch_size,
					size_t	pixel_start,
					size_t	pixel_end,
					T		*img_in,
					size_t	image_size,
					int		xdim,
					int		ydim,
					int		zdim,
					int		xshift,
					int		yshift,
					int		zshift);
//----------------------------------------------------------------------------
/*void probRatio( int       blockIdx_x, 
				int       threadIdx_x,
				XFLOAT *d_Mccf,
				XFLOAT *d_Mpsi,
				XFLOAT *d_Maux,
				XFLOAT *d_Mmean,
				XFLOAT *d_Mstddev,
				size_t image_size,
				XFLOAT  normfft,
				XFLOAT  sum_ref_under_circ_mask,
				XFLOAT sum_ref2_under_circ_mask,
				XFLOAT expected_Pratio,
				int NpsiThisBatch,
				int startPsi,
				int totalPsis);

void rotateOnly(int              blockIdx_x, 
				int              blockIdx_y, 
				int              threadIdx_x,
				ACCCOMPLEX     *d_Faux,
				XFLOAT           psi,
				AccProjectorKernel &projector,
	
void rotateAndCtf(  int              blockIdx_x, 
					int              blockIdx_y, 
					int              threadIdx_x,
					ACCCOMPLEX     *d_Faux,
					XFLOAT          *d_ctf,
					XFLOAT           psi,
					AccProjectorKernel &projector,
					int              startPsi = 0);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A
 *|
void convol_A(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A, writes to C
 *|
void convol_A(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				ACCCOMPLEX *d_C,
				size_t       image_size);

|*
 * Multiplies many complex arrays A (in-place) by a single B, pixel-by-pixel, after conjugating A
 *|
void batch_convol_A(int           blockIdx_x,
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					size_t        image_size);

|*
* Multiplies many complex arrays A (not in-place) by a single B, pixel-by-pixel, after conjugating A
*|
void batch_convol_A(int          blockIdx_x,
					int          threadIdx_x,
					ACCCOMPLEX *d_A,
					ACCCOMPLEX *d_B,
					ACCCOMPLEX *d_C,
					size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B
 *|
void convol_B(  int          blockIdx_x, 
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B, writes to C
 *|
void convol_B(  int       blockIdx_x, 
				int       threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				ACCCOMPLEX *d_C,
				size_t    image_size);
|*
 * Multiplies many complex arrays A (in-place) by a single one B, pixel-by-pixel, after conjugating B
 *|
void batch_convol_B(int           blockIdx_x, 
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					size_t        image_size);
|*
 * Multiplies scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 *|
template <typename T>
void cpu_kernel_multi( T   *A,
			T   *OUT,
			T    S,
			size_t   image_size);
*/
/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
template <typename T>
void cpu_kernel_multi( T   *A,
			T    S,
			size_t     image_size);
/*
 * Multiplies scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
template <typename T>
void cpu_kernel_multi( T *A,
			T *B,
			T *OUT,
			T  S,
			size_t   image_size);
/*
void finalizeMstddev(   int       blockIdx_x, 
						int       threadIdx_x,
						XFLOAT   *Mstddev,
						XFLOAT   *aux,
						XFLOAT    S,
						size_t       image_size);

|*
 * In place squares array in place
 *
 *  A[i] = A[i]*A[i]
 *|
void square(int       blockIdx_x, 
			int       threadIdx_x,
			XFLOAT   *A,
			size_t       image_size);
*/
/*
 * Casts on device so we can copy_to_host directly into a multidimarray.
 *
template <typename T1, typename T2 >
void cast(  int blockIdx_x,
			int threadIdx_x,
			T1 *IN,
			T2 *OUT,
			size_t size)
{
	size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
	if(pixel<size)
		OUT[pixel] = IN[pixel];
}
*/
template<bool do_highpass>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void kernel_frequencyPass( int grid_size, int block_size,
					ACCCOMPLEX *A,
					long int     ori_size,
					size_t       Xdim,
					size_t       Ydim,
					size_t       Zdim,
					XFLOAT       edge_low,
					XFLOAT       edge_width,
					XFLOAT       edge_high,
					XFLOAT       angpix,
					size_t       image_size)
{
#ifdef DEBUG_CUDA
	if((size_t)grid_size*(size_t)block_size > (size_t)std::numeric_limits<int>::max())
		CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  grid_size*(size_t)block_size > (size_t)std::numeric_limits<int>::max()");
	if (image_size < 0)
		CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  image_size < 0");
#endif
	// TODO - why not a single loop over image_size pixels?
	for(int blk=0; blk<grid_size; blk++) {
		for(int tid=0; tid<block_size; tid++) {
			size_t texel = (size_t)tid + (size_t)blk*(size_t)block_size;

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
		} // tid
	} // blk
}

template<bool DATA3D>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
inline
#endif
void powerClass(int          gridSize,
				ACCCOMPLEX   *g_image,
				XFLOAT       *g_spectrum,
				size_t       image_size,
				size_t       spectrum_size,
				int          xdim,
				int          ydim,
				int          zdim,
				int          res_limit,
				XFLOAT      *g_highres_Xi2)
{
#ifdef DEBUG_CUDA
	if((size_t)gridSize*(size_t)POWERCLASS_BLOCK_SIZE > (size_t)std::numeric_limits<int>::max())
		CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  gridSize*(size_t)POWERCLASS_BLOCK_SIZE > (size_t)std::numeric_limits<int>::max()");
	if (image_size < 0)
		CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  image_size < 0");
#endif
	for(int bid=0; bid<gridSize; bid++)
	{
		XFLOAT normFaux;

		int x,y,xy,d;
		int xydim = xdim*ydim;
		bool coords_in_range(true);

		XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];
		for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++)
			s_highres_Xi2[tid] = (XFLOAT)0.;

		for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++){
			size_t voxel=(size_t)tid + (size_t)bid*(size_t)POWERCLASS_BLOCK_SIZE;
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
				size_t ires = (size_t)(sqrt((XFLOAT)d) + 0.5);
		#else
				size_t ires = (size_t)(sqrtf((XFLOAT)d) + 0.5f);
		#endif
				if((ires<spectrum_size) && coords_in_range)
				{
					normFaux = g_image[voxel].x*g_image[voxel].x + g_image[voxel].y*g_image[voxel].y;
					g_spectrum[ires] += normFaux;
					if(ires>=res_limit)
						s_highres_Xi2[tid] += normFaux;
				}
			}
		}

		for(int tid=1; tid<POWERCLASS_BLOCK_SIZE; tid++)
			s_highres_Xi2[0] += s_highres_Xi2[tid];

		g_highres_Xi2[0] += s_highres_Xi2[0];
	}
}

#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline void translatePixel(
		int x,
		int y,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT &real,
		XFLOAT &imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT s, c;
#ifdef ACC_DOUBLE_PRECISION
	sincos( x * tx + y * ty , &s, &c );
#else
	sincosf( x * tx + y * ty , &s, &c );
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline void translatePixel(
		int x,
		int y,
		int z,
		XFLOAT tx,
		XFLOAT ty,
		XFLOAT tz,
		XFLOAT &real,
		XFLOAT &imag,
		XFLOAT &tReal,
		XFLOAT &tImag)
{
	XFLOAT s, c;
#ifdef ACC_DOUBLE_PRECISION
	sincos( x * tx + y * ty + z * tz, &s, &c );
#else
	sincosf( x * tx + y * ty + z * tz, &s, &c );
#endif

	tReal = c * real - s * imag;
	tImag = c * imag + s * real;
}

// sincos lookup table optimization. Function translatePixel calls 
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions. 
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and 
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to 
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty). 
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline void  computeSincosLookupTable2D(unsigned long  trans_num,
                                        XFLOAT  *trans_x,
										XFLOAT  *trans_y,										
										int      xSize,
										int      ySize,
										XFLOAT  *sin_x,
										XFLOAT  *cos_x,
	                                    XFLOAT  *sin_y,
	                                    XFLOAT  *cos_y)
{
	for(unsigned long i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];

		for(int x=0; x<xSize; x++) {
		   unsigned long index = i * xSize + x;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( x * tx, &sin_x[index], &cos_x[index] );
#else
			sincosf( x * tx, &sin_x[index], &cos_x[index] );
#endif            
		}
		
		for(int y=0; y<ySize; y++) {
            unsigned long index = i * ySize + y;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( y * ty, &sin_y[index], &cos_y[index] );
#else
			sincosf( y * ty, &sin_y[index], &cos_y[index] );
#endif            
		}        
	}
}	                                    
				
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline void  computeSincosLookupTable3D(unsigned long  trans_num,
                                        XFLOAT  *trans_x,
										XFLOAT  *trans_y,
										XFLOAT  *trans_z,										
										int      xSize,
										int      ySize,
										int      zSize,
										XFLOAT  *sin_x,
										XFLOAT  *cos_x,
	                                    XFLOAT  *sin_y,
	                                    XFLOAT  *cos_y,
	                                    XFLOAT  *sin_z,
	                                    XFLOAT  *cos_z)
{	                                    
	for(unsigned long i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];
		XFLOAT tz = trans_z[i];

		for(int x=0; x<xSize; x++) {
		   unsigned long index = i * xSize + x;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( x * tx, &sin_x[index], &cos_x[index] );
#else
			sincosf( x * tx, &sin_x[index], &cos_x[index] );
#endif            
		}
		
		for(int y=0; y<ySize; y++) {
            unsigned long index = i * ySize + y;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( y * ty, &sin_y[index], &cos_y[index] );
#else
			sincosf( y * ty, &sin_y[index], &cos_y[index] );
#endif            
		}           
		
		for(int z=0; z<zSize; z++) {
			unsigned long index = i * zSize + z;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( z * tz, &sin_z[index], &cos_z[index] );
#else
			sincosf( z * tz, &sin_z[index], &cos_z[index] );
#endif            
		}        		
	}
}

template<bool invert>
void cpu_kernel_make_eulers_2D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *eulers,
		unsigned long orientation_num);

template<bool invert,bool doL, bool doR>
void cpu_kernel_make_eulers_3D(int grid_size, int block_size,
		XFLOAT *alphas,
		XFLOAT *betas,
		XFLOAT *gammas,
		XFLOAT *eulers,
		unsigned long orientation_num,
		XFLOAT *L,
		XFLOAT *R);

} // end of namespace CpuKernels

#endif /* HELPER_KERNELS_H_ */
