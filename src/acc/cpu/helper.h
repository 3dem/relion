#ifndef HELPER_KERNELS_H_
#define HELPER_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include "src/macros.h"
#include "src/acc/cpu/settings.h"
#include "src/acc/cpu/utilities.h"
#include "src/acc/cpu/projector.h"

namespace CpuKernels
{

#ifdef ACC_DOUBLE_PRECISION
#define FAILSAFE_PRIOR_MIN_LIM 1e-300
#else
#define FAILSAFE_PRIOR_MIN_LIM 1e-30
#endif

template<bool failsafe,typename weights_t>
void exponentiate_weights_coarse(
        int blockIdx_x,
        int blockIdx_y,
        int threadIdx_x,
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		weights_t *g_Mweight,
		XFLOAT avg_diff2,
		XFLOAT min_diff2,
        int nr_coarse_orient,
        int nr_coarse_trans)
{
	// blockid
	int bid  = blockIdx_x;
	int cid  = blockIdx_y;
	//threadid
	int tid = threadIdx_x;

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
void collect2jobs(  int     blockIdx_x,
					XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
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
	int bid = blockIdx_x;

	XFLOAT s_o_weights[SUMW_BLOCK_SIZE];
	XFLOAT s_thr_wsum_sigma2_offset[SUMW_BLOCK_SIZE];;
	XFLOAT s_thr_wsum_prior_offsetx_class[SUMW_BLOCK_SIZE];
	XFLOAT s_thr_wsum_prior_offsety_class[SUMW_BLOCK_SIZE];
	XFLOAT s_thr_wsum_prior_offsetz_class[SUMW_BLOCK_SIZE];
	
	long int pos = d_job_idx[bid];
    int job_size = d_job_num[bid];

    int pass_num = ceilfracf(job_size,SUMW_BLOCK_SIZE);
    
    for(int tid=0; tid<SUMW_BLOCK_SIZE; tid++) {
       	s_o_weights[tid]                    	= (XFLOAT)0.0;
       	s_thr_wsum_sigma2_offset[tid]       	= (XFLOAT)0.0;
       	s_thr_wsum_prior_offsetx_class[tid] 	= (XFLOAT)0.0;
       	s_thr_wsum_prior_offsety_class[tid] 	= (XFLOAT)0.0;
       	if(DATA3D)
      		s_thr_wsum_prior_offsety_class[tid] = (XFLOAT)0.0;
    }
    
        
    for (int pass = 0; pass < pass_num; pass++, pos+=SUMW_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
    {
        for(int tid=0; tid<SUMW_BLOCK_SIZE; tid++) {
        	if ((pass*SUMW_BLOCK_SIZE+tid)<job_size) // if there is a translation that needs to be done still for this thread
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
    
   	for(int tid=1; tid<SUMW_BLOCK_SIZE; tid++)
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
}

void exponentiate_weights_fine(
                   int blockIdx_x,
                   int threadIdx_x,    
                   XFLOAT *g_pdf_orientation,
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

void softMaskBackgroundValue(	int     blockIdx_x,
								int     threadIdx_x,
								int     gridDim_x,
								XFLOAT *vol,
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

void cosineFilter(	int      blockIdx_x,
					int      threadIdx_x,
					int      gridDim_x,
					XFLOAT *vol,
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

void translate2D(   int      blockIdx_x,
					int      threadIdx_x,
					XFLOAT * g_image_in,
					XFLOAT * g_image_out,
					int      image_size,
					int      xdim,
					int      ydim, //not used
					int      dx,
					int      dy);

void translate3D(   int      blockIdx_x,
					int      threadIdx_x,
					XFLOAT * g_image_in,
					XFLOAT * g_image_out,
					int      image_size,
					int      xdim,
					int      ydim,
					int      zdim, //not used
					int      dx,
					int      dy,
					int      dz);

//----------------------------------------------------------------------------
void centerFFT_2D(  int       blockIdx_x,
					int       blockIdx_y,
					int       threadIdx_x,
					XFLOAT   *img_in,
					int       image_size,
					int       xdim,
					int       ydim,
					int       xshift,
					int       yshift);

void centerFFT_3D(  int       blockIdx_x,
					int       blockIdx_y,
					int       threadIdx_x,
					XFLOAT   *img_in,
					int       image_size,
					int       xdim,
					int       ydim,
					int       zdim,
					int       xshift,
					int       yshift,
					int       zshift);
//----------------------------------------------------------------------------
void probRatio( int       blockIdx_x, 
				int       threadIdx_x,
				XFLOAT *d_Mccf,
				XFLOAT *d_Mpsi,
				XFLOAT *d_Maux,
				XFLOAT *d_Mmean,
				XFLOAT *d_Mstddev,
				int     image_size,
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
				ProjectorKernel &projector,
				int              startPsi);

void rotateAndCtf(  int              blockIdx_x, 
					int              blockIdx_y, 
					int              threadIdx_x,
					ACCCOMPLEX     *d_Faux,
					XFLOAT          *d_ctf,
					XFLOAT           psi,
					ProjectorKernel &projector,
					int              startPsi = 0);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A
 */
void convol_A(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				int          image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A, writes to C
 */
void convol_A(  int          blockIdx_x,
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				ACCCOMPLEX *d_C,
				int          image_size);

/*
 * Multiplies many complex arrays A (in-place) by a single B, pixel-by-pixel, after conjugating A
 */
void batch_convol_A(int           blockIdx_x,
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					int           image_size);

/*
* Multiplies many complex arrays A (not in-place) by a single B, pixel-by-pixel, after conjugating A
*/
void batch_convol_A(int          blockIdx_x,
					int          threadIdx_x,
					ACCCOMPLEX *d_A,
					ACCCOMPLEX *d_B,
					ACCCOMPLEX *d_C,
					int          image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B
 */
void convol_B(  int          blockIdx_x, 
				int          threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				int          image_size);

/*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B, writes to C
 */
void convol_B(  int       blockIdx_x, 
				int       threadIdx_x,
				ACCCOMPLEX *d_A,
				ACCCOMPLEX *d_B,
				ACCCOMPLEX *d_C,
				int      image_size);
/*
 * Multiplies many complex arrays A (in-place) by a single one B, pixel-by-pixel, after conjugating B
 */
void batch_convol_B(int           blockIdx_x, 
					int           threadIdx_x,
					ACCCOMPLEX  *d_A,
					ACCCOMPLEX  *d_B,
					int           image_size);
/*
 * Multiplies scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 */
void multi( int       blockIdx_x, 
			int       threadIdx_x,
			XFLOAT   *A,
			XFLOAT   *OUT,
			XFLOAT    S,
			int       image_size);

/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
void multi( int       blockIdx_x, 
			int       threadIdx_x,
			XFLOAT   *A,
			XFLOAT    S,
			int       image_size);
/*
 * Multiplies scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
void multi( int     blockIdx_x, 
			int     threadIdx_x,
			XFLOAT *A,
			XFLOAT *B,
			XFLOAT *OUT,
			XFLOAT  S,
			int     image_size);

void finalizeMstddev(   int       blockIdx_x, 
						int       threadIdx_x,
						XFLOAT   *Mstddev,
						XFLOAT   *aux,
						XFLOAT    S,
						int       image_size);

/*
 * In place squares array in place
 *
 *  A[i] = A[i]*A[i]
 */
void square(int       blockIdx_x, 
			int       threadIdx_x,
			XFLOAT   *A,
			int       image_size);

/*
 * Casts on device so we can copy_to_host directly into a multidimarray.
 */
template <typename T1, typename T2 >
void cast(  int blockIdx_x,
			int threadIdx_x,
			T1 *IN,
			T2 *OUT,
			int size)
{
	int pixel = threadIdx_x + blockIdx_x*BLOCK_SIZE;
	if(pixel<size)
		OUT[pixel] = IN[pixel];
}

template<bool do_highpass>
void frequencyPass( int          blockIdx_x,
					int          threadIdx_x,
					ACCCOMPLEX *A,
					long int     ori_size,
					size_t       Xdim,
					size_t       Ydim,
					size_t       Zdim,
					XFLOAT       edge_low,
					XFLOAT       edge_width,
					XFLOAT       edge_high,
					XFLOAT       angpix,
					int          image_size)
{
	int texel = threadIdx_x + blockIdx_x*BLOCK_SIZE;

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
void powerClass(int           blockIdx_x,
				ACCCOMPLEX *g_image,
				XFLOAT      *g_spectrum,
				int          image_size,
				int          spectrum_size,
				int          xdim,
				int          ydim,
				int          zdim,
				int          res_limit,
				XFLOAT      *g_highres_Xi2)
{
	int bid =  blockIdx_x;

	XFLOAT normFaux;

	int x,y,xy,d;
	int xydim = xdim*ydim;
	bool coords_in_range(true);
	
	XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];
    for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++)
        s_highres_Xi2[tid] = (XFLOAT)0.;
        
    for(int tid=0; tid<POWERCLASS_BLOCK_SIZE; tid++){
        int voxel=tid + bid*POWERCLASS_BLOCK_SIZE;
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
	    	int ires = (int)(sqrt((XFLOAT)d) + 0.5);
#else
	    	int ires = (int)(sqrtf((XFLOAT)d) + 0.5f);
#endif
	    	if((ires>0.f) && (ires<spectrum_size) && coords_in_range)
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
inline void  computeSincosLookupTable2D(int      trans_num,
                                        XFLOAT  *trans_x,
										XFLOAT  *trans_y,										
										int      xSize,
										int      ySize,
										XFLOAT  *sin_x,
										XFLOAT  *cos_x,
	                                    XFLOAT  *sin_y,
	                                    XFLOAT  *cos_y)
{
	for(int i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];

		for(int x=0; x<xSize; x++) {
		   int index = i * xSize + x;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( x * tx, &sin_x[index], &cos_x[index] );
#else
			sincosf( x * tx, &sin_x[index], &cos_x[index] );
#endif            
		}
		
		for(int y=0; y<ySize; y++) {
            int index = i * ySize + y;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( y * ty, &sin_y[index], &cos_y[index] );
#else
			sincosf( y * ty, &sin_y[index], &cos_y[index] );
#endif            
		}        
	}
}	                                    
				
inline void  computeSincosLookupTable3D(int      trans_num,
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
	for(int i=0; i<trans_num; i++) {
		XFLOAT tx = trans_x[i];
		XFLOAT ty = trans_y[i];
		XFLOAT tz = trans_z[i];

		for(int x=0; x<xSize; x++) {
		   int index = i * xSize + x;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( x * tx, &sin_x[index], &cos_x[index] );
#else
			sincosf( x * tx, &sin_x[index], &cos_x[index] );
#endif            
		}
		
		for(int y=0; y<ySize; y++) {
            int index = i * ySize + y;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( y * ty, &sin_y[index], &cos_y[index] );
#else
			sincosf( y * ty, &sin_y[index], &cos_y[index] );
#endif            
		}           
		
		for(int z=0; z<zSize; z++) {
			int index = i * zSize + z;
#ifdef ACC_DOUBLE_PRECISION
			sincos ( z * tz, &sin_z[index], &cos_z[index] );
#else
			sincosf( z * tz, &sin_z[index], &cos_z[index] );
#endif            
		}        		
	}
}
									
} // end of namespace CpuKernels

#endif /* HELPER_KERNELS_H_ */
