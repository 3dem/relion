#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include <vector>
#include <iostream>

__global__ void cuda_kernel_sumweightCoarse(  FLOAT *g_pdf_orientation,
									     	  FLOAT *g_pdf_offset,
									     	  FLOAT *g_Mweight,
									     	  FLOAT *g_thisparticle_sumweight,
									     	  FLOAT min_diff2,
									     	  int oversamples_orient,
									     	  int oversamples_trans,
									     	  int coarse_trans)
{
	__shared__ FLOAT s_sumweight[SUM_BLOCK_SIZE];
	// blockid
	int bid  = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	s_sumweight[tid]=0.;
	int c_iorient, f_iorient, c_itrans, f_itrans, pos, iorient = bid*SUM_BLOCK_SIZE+tid;

	// Bacause the partion of work is so arbitrarily divided in this kernel,
	// we need to do some brute idex work to get the correct indices.
	c_iorient = (iorient - (iorient % oversamples_orient)) / oversamples_orient; //floor(x/y) == (x-(x%y))/y  but less sensitive to x>>y and finite precision
	f_iorient = iorient % oversamples_orient;
	for (int itrans=0; itrans<(coarse_trans*oversamples_trans); itrans++)
	{
		c_itrans = ( itrans - (itrans % oversamples_trans))/ oversamples_trans; //floor(x/y) == (x-(x%y))/y  but less sensitive to x>>y and finite precision
		f_itrans = itrans % oversamples_trans;

		pos = c_iorient * (coarse_trans *oversamples_orient*oversamples_trans );
		pos += c_itrans * (oversamples_orient*oversamples_trans);
		pos += f_iorient*( oversamples_trans ) + f_itrans;

		if( g_Mweight[pos] < (FLOAT)0.0 ) //TODO Might be slow (divergent threads)
		{
			g_Mweight[pos] = (FLOAT)0.0;
		}
		else
		{
			FLOAT weight = g_pdf_orientation[c_iorient] * g_pdf_offset[c_itrans];          	// Same      for all threads - TODO: should be done once for all trans through warp-parallel execution
			FLOAT diff2 = g_Mweight[pos] - min_diff2;								// Different for all threads
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
				// TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn

			// Store the weight
			g_Mweight[pos] = weight; // TODO put in shared mem

			// Reduce weights for sum of all weights
			s_sumweight[tid] += weight;
		}
	}

	__syncthreads();
	// Further reduction of all samples in this block
	for(int j=(SUM_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
		{
			s_sumweight[tid] += s_sumweight[tid+j];
		}
	}
	__syncthreads();
	g_thisparticle_sumweight[bid]=s_sumweight[0];
}


/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
__global__ void cuda_kernel_sumweightFine(    FLOAT *g_pdf_orientation,
									     	  FLOAT *g_pdf_offset,
									     	  FLOAT *g_weights,
									     	  FLOAT *g_thisparticle_sumweight,
									     	  FLOAT min_diff2,
									     	  int oversamples_orient,
									     	  int oversamples_trans,
									     	  unsigned long *d_rot_id,
									     	  unsigned long *d_trans_idx,
									     	  unsigned long *d_job_idx,
									     	  unsigned long *d_job_num,
									     	  long int job_num)
{
	__shared__ FLOAT s_sumweight[SUM_BLOCK_SIZE];
	__shared__ FLOAT s_weights[SUM_BLOCK_SIZE];

	// blockid
	int bid  = blockIdx.x;
	//threadid
	int tid = threadIdx.x;

	s_sumweight[tid]=0.;

	long int jobid = bid*SUM_BLOCK_SIZE+tid;

	if (jobid<job_num)
	{
		long int pos = d_job_idx[jobid];
		// index of comparison
		long int ix =  d_rot_id[   pos];   // each thread gets its own orient...
		long int iy = d_trans_idx[ pos];   // ...and it's starting trans...
		long int in =  d_job_num[jobid];    // ...AND the number of translations to go through

		int c_iorient, f_iorient, c_itrans, f_itrans, iorient = bid*SUM_BLOCK_SIZE+tid;

		// Bacause the partion of work is so arbitrarily divided in this kernel,
		// we need to do some brute idex work to get the correct indices.
		for (int itrans=0; itrans < in; itrans++, iy++)
		{
			c_itrans = ( iy - (iy % oversamples_trans))/ oversamples_trans; //floor(x/y) == (x-(x%y))/y  but less sensitive to x>>y and finite precision
			f_itrans = iy % oversamples_trans;

			FLOAT prior = g_pdf_orientation[ix] * g_pdf_offset[c_itrans];          	// Same      for all threads - TODO: should be done once for all trans through warp-parallel execution
			FLOAT diff2 = g_weights[pos+itrans] - min_diff2;								// Different for all threads
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

			// Reduce weights for sum of all weights
			s_sumweight[tid] += s_weights[tid];
		}
	}
	else
	{
		s_sumweight[tid]=0.;
	}

	__syncthreads();
	// Further reduction of all samples in this block
	for(int j=(SUM_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
			s_sumweight[tid] += s_sumweight[tid+j];
	}
	__syncthreads();
	g_thisparticle_sumweight[bid]=s_sumweight[0];
}

__global__ void cuda_kernel_collect2(	FLOAT *g_oo_otrans_x,          // otrans-size -> make const
										FLOAT *g_oo_otrans_y,          // otrans-size -> make const
										FLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
										FLOAT *g_Mweight,
										FLOAT op_significant_weight,    // TODO Put in const
										FLOAT op_sum_weight,            // TODO Put in const
										int   coarse_trans,
										int   oversamples_trans,
										int   oversamples_orient,
										int   oversamples,
										bool  do_ignore_pdf_direction,
										FLOAT *g_weights,
										FLOAT *g_thr_wsum_prior_offsetx_class,
										FLOAT *g_thr_wsum_prior_offsety_class,
										FLOAT *g_thr_wsum_sigma2_offset
										)
{
	// objects reduced in this kernel, which need to be further reduced for all blocks
	// after the kernel has finished. Notice that all reductions span fine sampling =>
	// a block can treat all fine samples in a coarse orientation and output a single
	// floating point value for each reduction. We do however list the dimension of
	// post-kernel reduction for all reductions here:
	__shared__ FLOAT                      s_weights[SUM_BLOCK_SIZE];
	__shared__ FLOAT s_thr_wsum_prior_offsetx_class[SUM_BLOCK_SIZE];
	__shared__ FLOAT s_thr_wsum_prior_offsety_class[SUM_BLOCK_SIZE];
	__shared__ FLOAT       s_thr_wsum_sigma2_offset[SUM_BLOCK_SIZE];

	int ex  = blockIdx.x * gridDim.y + blockIdx.y;            // coarse orientation
	int tid = threadIdx.x;
	int pass_num = ceil((float)oversamples / (float)SUM_BLOCK_SIZE);
	//Where to start in g_Mweight to find all data for this *coarse* orientation
	long int ref_Mweight_idx = ex * ( coarse_trans*oversamples );

	int iover_trans = tid % oversamples_trans;
	int iover_rot = floor((float) tid / (float)oversamples_trans);
	s_weights[tid]                      = 0.0f;
	s_thr_wsum_prior_offsetx_class[tid] = 0.0f;
	s_thr_wsum_prior_offsety_class[tid] = 0.0f;
	s_thr_wsum_sigma2_offset[tid]       = 0.0f;

	// Go over all (21 typically) *coarse* translations, reducing in place
	for (int itrans=0; itrans<coarse_trans; itrans++)
	{
		//Where to start in g_Mweights to find all fine samples for this *coarse* translation
		int pos = ref_Mweight_idx + itrans*oversamples + iover_rot*oversamples_trans + iover_trans;
		for (int pass = 0; pass < pass_num; pass++, pos+=SUM_BLOCK_SIZE)
		{
			FLOAT weight = g_Mweight[pos];
			if( weight >= op_significant_weight ) //TODO Might be slow (divergent threads)
				weight /= op_sum_weight;
			else
				weight = 0.0f;

			s_weights[tid] += weight;
			s_thr_wsum_prior_offsetx_class[tid] +=    weight * g_oo_otrans_x[iover_trans + itrans*oversamples_trans];    // precalc otrans_y, only overtrans-size => const
			s_thr_wsum_prior_offsety_class[tid] +=    weight * g_oo_otrans_y[iover_trans + itrans*oversamples_trans];    // precalc otrans_y, only overtrans-size => const
			s_thr_wsum_sigma2_offset[tid] += weight * g_myp_oo_otrans_x2y2z2[iover_trans + itrans*oversamples_trans];    // precalc x2y2z2,   only overtrans-size => const
		}
	}
	// Reduction of all fine samples in this coarse orientation
	for(int j=(SUM_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
		{
			s_weights[tid]                      += s_weights[tid+j];
			s_thr_wsum_prior_offsetx_class[tid] += s_thr_wsum_prior_offsetx_class[tid+j];
			s_thr_wsum_prior_offsety_class[tid] += s_thr_wsum_prior_offsety_class[tid+j];
			s_thr_wsum_sigma2_offset[tid]       += s_thr_wsum_sigma2_offset[tid+j];
		}
		__syncthreads();
	}
	// write pre-reduced (for all fine samples and itrans) to global mem.
	g_weights[ex]			           = s_weights[0];
	g_thr_wsum_prior_offsetx_class[ex] = s_thr_wsum_prior_offsetx_class[0];
	g_thr_wsum_prior_offsety_class[ex] = s_thr_wsum_prior_offsety_class[0];
	g_thr_wsum_sigma2_offset[ex]       = s_thr_wsum_sigma2_offset[0];
}

__global__ void cuda_kernel_collect2jobs(	FLOAT *g_oo_otrans_x,          // otrans-size -> make const
										FLOAT *g_oo_otrans_y,          // otrans-size -> make const
										FLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
										FLOAT *g_i_weights,
										FLOAT op_significant_weight,    // TODO Put in const
										FLOAT op_sum_weight,            // TODO Put in const
										int   coarse_trans,
										int   oversamples_trans,
										int   oversamples_orient,
										int   oversamples,
										bool  do_ignore_pdf_direction,
										FLOAT *g_o_weights,
										FLOAT *g_thr_wsum_prior_offsetx_class,
										FLOAT *g_thr_wsum_prior_offsety_class,
										FLOAT *g_thr_wsum_sigma2_offset,
								     	unsigned long *d_rot_idx,
								     	unsigned long *d_trans_idx,
								     	unsigned long *d_job_idx,
								     	unsigned long *d_job_num
								     	)
{
	// blockid
	int bid  =blockIdx.x * gridDim.y + blockIdx.y;
	//threadid
	int tid = threadIdx.x;

	__shared__ FLOAT                    s_o_weights[SUM_BLOCK_SIZE];
	__shared__ FLOAT s_thr_wsum_prior_offsetx_class[SUM_BLOCK_SIZE];
	__shared__ FLOAT s_thr_wsum_prior_offsety_class[SUM_BLOCK_SIZE];
	__shared__ FLOAT       s_thr_wsum_sigma2_offset[SUM_BLOCK_SIZE];
	s_o_weights[tid]                    = (FLOAT)0.0;
	s_thr_wsum_prior_offsetx_class[tid] = (FLOAT)0.0;
	s_thr_wsum_prior_offsety_class[tid] = (FLOAT)0.0;
	s_thr_wsum_sigma2_offset[tid]       = (FLOAT)0.0;

	long int pos = d_job_idx[bid];
    int job_size = d_job_num[bid];
	pos += tid;	   					// pos is updated to be thread-resolved

    int pass_num = ceil((float)job_size / (float)SUM_BLOCK_SIZE);
    for (int pass = 0; pass < pass_num; pass++, pos+=SUM_BLOCK_SIZE) // loop the available warps enough to complete all translations for this orientation
    {
    	if ((pass*SUM_BLOCK_SIZE+tid)<job_size) // if there is a translation that needs to be done still for this thread
    	{
			// index of comparison
			long int iy = d_trans_idx[pos];              // ...and its own trans...

			FLOAT weight = g_i_weights[pos];
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
	for(int j=(SUM_BLOCK_SIZE/2); j>0; j/=2)
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


__global__ void cuda_kernel_reduce_wdiff2s(FLOAT *g_wdiff2s_parts,
										   long int orientation_num,
										   int image_size,
										   int current_block_num)
{
	unsigned long bid = blockIdx.y*gridDim.x + blockIdx.x;
	unsigned tid = threadIdx.x;
	unsigned pass_num(ceilf((float)image_size/(float)BLOCK_SIZE)),pixel;
	if((current_block_num+bid)<orientation_num)
	{
		for (unsigned pass = 0; pass < pass_num; pass++)
		{
			pixel = pass * BLOCK_SIZE + tid;
			if(pixel<image_size)
				g_wdiff2s_parts[bid*image_size+pixel] += g_wdiff2s_parts[(current_block_num+bid)*image_size+pixel];
		}
	}
}
