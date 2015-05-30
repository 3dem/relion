#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include <vector>
#include <iostream>

__global__ void cuda_kernel_sumweight_oversampling(	FLOAT *g_pdf_orientation,
													FLOAT *g_pdf_offset,
													FLOAT *g_Mweight,
													FLOAT *g_thisparticle_sumweight,
													FLOAT min_diff2,
													int translation_num,
													int oversamples)
{
	__shared__ FLOAT s_sumweight[SUM_BLOCK_SIZE];
	// blockid
	int ex  = blockIdx.x * gridDim.y + blockIdx.y;
	//threadid
	int tid = threadIdx.x;
	s_sumweight[tid]=0;

	// passes to take care of all fine samples in a coarse sample
	int pass_num = ceil((float)oversamples / (float)SUM_BLOCK_SIZE);
	//Where to start in g_Mweight to find all data for this *coarse* orientation
	long int ref_Mweight_idx = ex * ( translation_num*oversamples );

	// Go over all *coarse* translations, reducing in place
	for (int itrans=0; itrans<translation_num; itrans++)
	{
		//Where to start in g_Mweights to find all fine samples for this *coarse* translation
		int pos = ref_Mweight_idx + itrans*oversamples + tid;
		for (int pass = 0; pass < pass_num; pass++, pos+=SUM_BLOCK_SIZE)
		{
			if( g_Mweight[pos] < 0.0f ) //TODO Might be slow (divergent threads)
			{
				g_Mweight[pos] = 0.0f;
			}
			else
			{
				FLOAT weight = g_pdf_orientation[ex] * g_pdf_offset[itrans];          	// Same      for all threads - TODO: should be done once for all trans through warp-parallel execution
				FLOAT diff2 = g_Mweight[pos] - min_diff2;								// Different for all threads
				// next line because of numerical precision of exp-function
				if (diff2 > 700.0f)
					weight = 0.0f;
				else weight *= exp(-diff2);  // TODO: use tabulated exp function? / Sjors  TODO: exp, expf, or __exp in CUDA? /Bjorn

				// Store the weight for each fine sample in this coarse pair
				g_Mweight[pos] = weight; // TODO put in shared mem

				// Reduce weights for each fine sample in this coarse pair
				s_sumweight[tid] += weight;
			}
		}
	}
	// Reduction of all fine samples in this coarse orientation
	for(int j=(SUM_BLOCK_SIZE/2); j>0; j/=2)
	{
		if(tid<j)
		{
			s_sumweight[tid] += s_sumweight[tid+j];
		}
	}
	__syncthreads();
	g_thisparticle_sumweight[ex]=s_sumweight[0];
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
