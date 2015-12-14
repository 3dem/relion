#include "src/gpu_utils/cuda_device_utils.cuh"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_settings.h"

__global__ void cuda_kernel_exponentiate_weights_coarse(
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		XFLOAT *g_Mweight,
		XFLOAT min_diff2,
		int nr_coarse_orient,
		int nr_coarse_trans)
{
	// blockid
	int bid  = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	int pos, iorient = bid*SUMW_BLOCK_SIZE+tid;

	XFLOAT weight;
	if(iorient<nr_coarse_orient)
	{
		for (int itrans=0; itrans<nr_coarse_trans; itrans++)
		{
			pos = iorient * nr_coarse_trans + itrans;
			XFLOAT diff2 = g_Mweight[pos] - min_diff2;
			if( diff2 < (XFLOAT)0.0 ) //TODO Might be slow (divergent threads)
				diff2 = (XFLOAT)0.0;
			else
			{
				weight = g_pdf_orientation[iorient] * g_pdf_offset[itrans];          	// Same for all threads - TODO: should be done once for all trans through warp-parallel execution

				// next line because of numerical precision of exp-function
#if defined(CUDA_DOUBLE_PRECISION)
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
				diff2=weight;
				// TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn
			}

			// Store the weight
			g_Mweight[pos] = diff2; // TODO put in shared mem
		}
	}
}


/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
__global__ void cuda_kernel_exponentiate_weights_fine(
		XFLOAT *g_pdf_orientation,
		XFLOAT *g_pdf_offset,
		XFLOAT *g_weights,
		XFLOAT min_diff2,
		int oversamples_orient,
		int oversamples_trans,
		unsigned long *d_rot_id,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num,
		long int job_num)
{
	__shared__ XFLOAT s_weights[SUMW_BLOCK_SIZE];

	// blockid
	int bid  = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	long int jobid = bid*SUMW_BLOCK_SIZE+tid;

	if (jobid<job_num)
	{
		long int pos = d_job_idx[jobid];
		// index of comparison
		long int ix =  d_rot_id[   pos];   // each thread gets its own orient...
		long int iy = d_trans_idx[ pos];   // ...and it's starting trans...
		long int in =  d_job_num[jobid];    // ...AND the number of translations to go through

		int c_itrans;//, iorient = bid*SUM_BLOCK_SIZE+tid; //, f_itrans;

		// Bacause the partion of work is so arbitrarily divided in this kernel,
		// we need to do some brute idex work to get the correct indices.
		for (int itrans=0; itrans < in; itrans++, iy++)
		{
			c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans; //floor(x/y) == (x-(x%y))/y  but less sensitive to x>>y and finite precision
//			f_itrans = iy % oversamples_trans;

			XFLOAT prior = g_pdf_orientation[ix] * g_pdf_offset[c_itrans];          	// Same      for all threads - TODO: should be done once for all trans through warp-parallel execution
			XFLOAT diff2 = g_weights[pos+itrans] - min_diff2;								// Different for all threads
			// next line because of numerical precision of exp-function
	#if defined(CUDA_DOUBLE_PRECISION)
				if (diff2 > 700.)
					s_weights[tid] = 0.;
				else
					s_weights[tid] = prior * exp(-diff2);
	#else
				if (diff2 > 88.)
					s_weights[tid] = 0.;
				else
					s_weights[tid] = prior * expf(-diff2);
	#endif
				// TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn
			// Store the weight
			g_weights[pos+itrans] = s_weights[tid]; // TODO put in shared mem
		}
	}
}

__global__ void cuda_kernel_collect2jobs(	XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
										XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
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

	__shared__ XFLOAT                    s_o_weights[SUMW_BLOCK_SIZE];
	__shared__ XFLOAT s_thr_wsum_prior_offsetx_class[SUMW_BLOCK_SIZE];
	__shared__ XFLOAT s_thr_wsum_prior_offsety_class[SUMW_BLOCK_SIZE];
	__shared__ XFLOAT       s_thr_wsum_sigma2_offset[SUMW_BLOCK_SIZE];
	s_o_weights[tid]                    = (XFLOAT)0.0;
	s_thr_wsum_prior_offsetx_class[tid] = (XFLOAT)0.0;
	s_thr_wsum_prior_offsety_class[tid] = (XFLOAT)0.0;
	s_thr_wsum_sigma2_offset[tid]       = (XFLOAT)0.0;

	long int pos = d_job_idx[bid];
    int job_size = d_job_num[bid];
	pos += tid;	   					// pos is updated to be thread-resolved

    int pass_num = ceilfracf(job_size,SUMW_BLOCK_SIZE);
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
				weight = 0.0f;

			s_o_weights[tid] += weight;
			s_thr_wsum_prior_offsetx_class[tid] += weight *          g_oo_otrans_x[iy];
			s_thr_wsum_prior_offsety_class[tid] += weight *          g_oo_otrans_y[iy];
			s_thr_wsum_sigma2_offset[tid]       += weight * g_myp_oo_otrans_x2y2z2[iy];
    	}
    }
    // Reduction of all treanslations this orientation
	for(int j=(SUMW_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
		{
			s_o_weights[tid]                    += s_o_weights[tid+j];
			s_thr_wsum_prior_offsetx_class[tid] += s_thr_wsum_prior_offsetx_class[tid+j];
			s_thr_wsum_prior_offsety_class[tid] += s_thr_wsum_prior_offsety_class[tid+j];
			s_thr_wsum_sigma2_offset[tid]       += s_thr_wsum_sigma2_offset[tid+j];
		}
		__syncthreads();
	}
	g_o_weights[bid]			           = s_o_weights[0];
	g_thr_wsum_prior_offsetx_class[bid] = s_thr_wsum_prior_offsetx_class[0];
	g_thr_wsum_prior_offsety_class[bid] = s_thr_wsum_prior_offsety_class[0];
	g_thr_wsum_sigma2_offset[bid]       = s_thr_wsum_sigma2_offset[0];
}
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
												XFLOAT cosine_width	)
{

		int tid = threadIdx.x;

//		vol.setXmippOrigin(); // sets xinit=xdim , also for y z
		XFLOAT r, raisedcos;

		__shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT    partial_sum[SOFTMASK_BLOCK_SIZE];
		__shared__ XFLOAT partial_sum_bg[SOFTMASK_BLOCK_SIZE];

		XFLOAT sum_bg_total = 0.f;

		long int texel_pass_num = ceilfracf(vol_size,SOFTMASK_BLOCK_SIZE);
		int texel = tid;

		partial_sum[tid]=0.f;
		partial_sum_bg[tid]=0.f;
		if (do_Mnoise)
		{
			for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
			{
				XFLOAT x,y,z;
				if(texel<vol_size)
				{
					img_pixels[tid]=__ldg(&vol[texel]);

					z = 0.f;// floor( (float) texel                  / (float)((xdim)*(ydim)));
					y = floor( (float)(texel-z*(xdim)*(ydim)) / (float) xdim );
					x = texel - z*(xdim)*(ydim) - y*xdim;

	//				z-=zinit;
					y-=yinit;
					x-=xinit;

					r = sqrt(x*x + y*y);// + z*z);

					if (r < radius)
						continue;
					else if (r > radius_p)
					{
						partial_sum[tid]    += 1.f;
						partial_sum_bg[tid] += img_pixels[tid];
					}
					else
					{
						raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
						partial_sum[tid] += raisedcos;
						partial_sum_bg[tid] += raisedcos * img_pixels[tid];
					}
				}
			}
		}

		__syncthreads();
		for(int j=(SOFTMASK_BLOCK_SIZE/2); j>0; j/=2)
		{
			if(tid<j)
			{
				partial_sum[tid] += partial_sum[tid+j];
				partial_sum_bg[tid] += partial_sum_bg[tid+j];
			}
			__syncthreads();
		}

		sum_bg_total  = partial_sum_bg[0] / partial_sum[0];


		texel = tid;
		for (int pass = 0; pass < texel_pass_num; pass++, texel+=SOFTMASK_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
		{
			XFLOAT x,y,z;
			if(texel<vol_size)
			{
				img_pixels[tid]=__ldg(&vol[texel]);

				z = 0.f;// floor( (float) texel                  / (float)((xdim)*(ydim)));
				y = floor( (float)(texel-z*(xdim)*(ydim)) / (float)  xdim         );
				x = texel - z*(xdim)*(ydim) - y*xdim;

//				z-=zinit;
				y-=yinit;
				x-=xinit;

				r = sqrt(x*x + y*y);// + z*z);

				if (r < radius)
					continue;
				else if (r > radius_p)
					img_pixels[tid]=sum_bg_total;
				else
				{
					raisedcos = 0.5f + 0.5f * cospif((radius_p - r) / cosine_width );
					img_pixels[tid]= img_pixels[tid]*(1-raisedcos) + sum_bg_total*raisedcos;
				}
				vol[texel]=img_pixels[tid];
			}

		}
}

__global__ void cuda_kernel_centerFFT_2D(XFLOAT *img_in,
										 XFLOAT *img_out,
										 long int image_size,
										 long int xdim,
										 long int ydim,
										 long int xshift,
										 long int yshift)
{
	long int pixel = threadIdx.x + blockIdx.x*CFTT_BLOCK_SIZE;
//	int pixel_pass_num = ceilfracf(image_size, CFTT_BLOCK_SIZE);

//	for (int pass = 0; pass < pixel_pass_num; pass++, pixel+=CFTT_BLOCK_SIZE)
//	{
		if(pixel<image_size)
		{
			int y = floorfracf(pixel,xdim);
			int x = pixel % xdim;				// also = pixel - y*xdim, but this depends on y having been calculated, i.e. serial evaluation

			int yp = y + yshift;
			if (yp < 0)
				yp += ydim;
			else if (yp >= ydim)
				yp -= ydim;

			int xp = x + xshift;
			if (xp < 0)
				xp += xdim;
			else if (xp >= xdim)
				xp -= xdim;

			long int n_pixel = yp*xdim + xp;

			img_out[n_pixel] = __ldg(&img_in[pixel]);
		}
//	}
}

__global__ void cuda_kernel_probRatio(  XFLOAT *d_Mccf,
										XFLOAT *d_Mpsi,
										XFLOAT *d_Maux,
										XFLOAT *d_Mmean,
										XFLOAT *d_Mstddev,
										long int image_size,
										XFLOAT normfft,
										XFLOAT sum_ref_under_circ_mask,
										XFLOAT sum_ref2_under_circ_mask,
										XFLOAT expected_Pratio,
										XFLOAT psi)
{
	/* PLAN TO:
	 *
	 * 1) Pre-filter
	 * 		d_Mstddev[i] = 1 / (2*d_Mstddev[i])   ( if d_Mstddev[pixel] > 1E-10 )
	 * 		d_Mstddev[i] = 1    				  ( else )
	 *
	 * 2) Set
	 * 		sum_ref2_under_circ_mask /= 2.
	 *
	 * 3) Total expression becomes
	 * 		diff2 = ( exp(k) - 1.f ) / (expected_Pratio - 1.f)
	 * 	  where
	 * 	  	k = (normfft * d_Maux[pixel] + d_Mmean[pixel] * sum_ref_under_circ_mask)*d_Mstddev[i] + sum_ref2_under_circ_mask
	 *
	 */

	long int pixel = threadIdx.x + blockIdx.x*PROBRATIO_BLOCK_SIZE;

	if(pixel<image_size)
	{
		XFLOAT diff2 = normfft * d_Maux[pixel];
		diff2 += d_Mmean[pixel] * sum_ref_under_circ_mask;

		if (d_Mstddev[pixel] > 1E-10)
			diff2 /= d_Mstddev[pixel];
		diff2 += sum_ref2_under_circ_mask;

#if defined(CUDA_DOUBLE_PRECISION)
		diff2 = exp(-diff2 / 2.); // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)
#else
		diff2 = expf(-diff2 / 2.f);
#endif

		// Store fraction of (1 - probability-ratio) wrt  (1 - expected Pratio)
		diff2 = (diff2 - 1.f) / (expected_Pratio - 1.f);
		if (diff2 > d_Mccf[pixel])
		{
			d_Mccf[pixel] = diff2;
			d_Mpsi[pixel] = psi;
		}
	}

}

//__global__ void cuda_kernel_rotateAndCtf( CUDACOMPLEX *d_Faux,
//						  	  	  	  	  XFLOAT *d_ctf,
//						  	  	  	  	  XFLOAT psi,
//						  	  			  CudaProjectorKernel projector)
//{
//
//	int image_size=projector.imgX*projector.imgY;
//	long int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
//	if(pixel<image_size)
//	{
//		int y = floorfracf(pixel,projector.imgX);
//		int x = pixel % projector.imgX;
//
//		if (y > projector.maxR)
//		{
//			if (y >= projector.imgY - projector.maxR)
//				y = y - projector.imgY;
//			else
//				x = projector.maxR;
//		}
//
//		XFLOAT sa, ca;
//		sincos(psi, &sa, &ca);
//		CUDACOMPLEX val;
//
//		projector.project2Dmodel(	 x,y,
//									 ca,
//									-sa,
//									 sa,
//									 ca,
//									 val.x,val.y);
//		d_Faux[pixel].x =val.x*d_ctf[pixel];
//		d_Faux[pixel].y =val.y*d_ctf[pixel];
//
//	}
//}

__global__ void cuda_kernel_convol(	 CUDACOMPLEX *d_A,
									 CUDACOMPLEX *d_B,
									 int image_size)
{
	long int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
	{
		XFLOAT tr =d_A[pixel].x;
		XFLOAT ti =d_A[pixel].y;
		d_A[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
		d_A[pixel].y = - ti*d_B[pixel].x + tr*d_B[pixel].y;
	}
}
