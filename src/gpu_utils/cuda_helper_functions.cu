#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <signal.h>

/*
 * This assisting function goes over the orientations determined as significant for this image, and checks
 * which translations should be included in the list of those which differences will be calculated for.
 *
 * Any contiguous translations with a shared orientation are grouped together into a "job" which is supplied
 * to the difference kernel. If there are more contiguous translations than the specified PROJDIFF_CHUNK_SIZE,
 * these are split into separate jobs, to increase paralllelism at the cost of redundant memory reads.
 */
long int makeJobsForDiff2Fine(
		OptimisationParamters &op,  SamplingParameters &sp,
		long int orientation_num, long int translation_num,
		ProjectionParams &FineProjectionData,
		std::vector< long unsigned > &iover_transes,
		std::vector< long unsigned > &ihiddens,
		long int nr_over_orient, long int nr_over_trans, int ipart,
		IndexedDataArray &FPW, // FPW=FinePassWeights
		IndexedDataArrayMask &dataMask)
{
	long int w_base = dataMask.firstPos, w(0), k(0);
	// be on the safe side with the jobArrays: make them as large as they could possibly be
	// (this will be reduced at exit of this function)
	dataMask.setNumberOfJobs(orientation_num*translation_num);
	dataMask.setNumberOfWeights(orientation_num*translation_num);
	dataMask.jobOrigin.host_alloc();
	dataMask.jobExtent.host_alloc();

	dataMask.jobOrigin[k]=0;
	for (long unsigned i = 0; i < orientation_num; i++)
	{
		dataMask.jobExtent[k]=0;
		int tk=0;
		long int iover_rot = FineProjectionData.iover_rots[i];
		for (long unsigned j = 0; j < translation_num; j++)
		{
			long int iover_trans = iover_transes[j];
			long int ihidden = FineProjectionData.iorientclasses[i] * sp.nr_trans + ihiddens[j];

			if(DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden)==1)
			{
				FPW.rot_id[w_base+w] = FineProjectionData.iorientclasses[i] % (sp.nr_dir*sp.nr_psi); 	// where to look for priors etc
				FPW.rot_idx[w_base+w] = i;					// which rot for this significant task
				FPW.trans_idx[w_base+w] = j;					// which trans       - || -
				FPW.ihidden_overs[w_base+w]= (ihidden * nr_over_orient + iover_rot) * nr_over_trans + iover_trans;

				if(tk>=PROJDIFF_CHUNK_SIZE)
				{
					tk=0;             // reset counter
					k++;              // use new element
					dataMask.jobOrigin[k]=w;
					dataMask.jobExtent[k]=0;   // prepare next element for ++ incrementing
				}
				tk++;                 		   // increment limit-checker
				dataMask.jobExtent[k]++;       // increment number of transes this job
				w++;
			}
			else if(tk!=0) 		  // start a new one with the same rotidx - we expect transes to be sequential.
			{
				tk=0;             // reset counter
				k++;              // use new element
				dataMask.jobOrigin[k]=w;
				dataMask.jobExtent[k]=0;   // prepare next element for ++ incrementing
			}
		}
		if(tk>0) // use new element (if tk==0) then we are currently on an element with no signif, so we should continue using this element
		{
			k++;
			dataMask.jobOrigin[k]=w;
			dataMask.jobExtent[k]=0;
		}
	}
	if(dataMask.jobExtent[k]!=0) // if we started putting somehting in last element, then the count is one higher than the index
		k+=1;

	dataMask.setNumberOfJobs(k);
	dataMask.setNumberOfWeights(w);
//	if(dataMask.weightNum>0)
//	{
//		dataMask.jobOrigin.device_alloc();
//		dataMask.jobExtent.device_alloc();
//	}
	return(w);
}

/*
 * This assisting function goes over the weight-array and groups all weights with shared
 * orientations into 'jobs' which are fed into the collect-kenrel, which reduces all translations
 * with computed differences into a reduced object to be back-projected.
 */
int  makeJobsForCollect(IndexedDataArray &FPW, IndexedDataArrayMask &dataMask, unsigned long NewJobNum) // FPW=FinePassWeights
{
	// reset the old (diff2Fine) job-definitions
//	dataMask.jobOrigin.free_host();
//    dataMask.jobOrigin.free_device();
//    dataMask.jobExtent.free_host();
//    dataMask.jobExtent.free_device();
    dataMask.setNumberOfJobs(NewJobNum);
//    dataMask.jobOrigin.host_alloc();
//    dataMask.jobExtent.host_alloc();

	long int jobid=0;
	dataMask.jobOrigin[jobid]=0;
	dataMask.jobExtent[jobid]=1;
	long int crot =FPW.rot_idx[jobid]; // set current rot
	for(long int n=1; n<FPW.rot_idx.size; n++)
	{
		if(FPW.rot_idx[n]==crot)
		{
			dataMask.jobExtent[jobid]++;
		}
		else
		{
			jobid++;
			dataMask.jobExtent[jobid]=1;
			dataMask.jobOrigin[jobid]=n;
			crot=FPW.rot_idx[n];
		}
	}
	dataMask.setNumberOfJobs(jobid+1); // because max index is one less than size
//	dataMask.jobOrigin.put_on_device();
//	dataMask.jobExtent.put_on_device();

	return (jobid+1);
}

/*
 * Splits a designated number of blocks to be run
 * in a cuda-kernel into two components of a dim3.
 */
dim3 splitCudaBlocks(long int block_num, bool doForceEven)
{
	unsigned int orient1, orient2;
	if(block_num>65535)
	{
		orient1 = ceil(sqrt(block_num));
		if(doForceEven)
			orient1 += 	(orient1 % 2);
		orient2 = orient1;
	}
	else
	{
		orient1 = block_num;
		orient2 = 1;
	}
	dim3 block_dim(orient1,orient2);
	return(block_dim);
}

/*
 * Maps weights to a decoupled indexing of translations and orientations
 */
void mapWeights(
		unsigned long orientation_start,
		XFLOAT *mapped_weights,
		unsigned orientation_num,
		unsigned long idxArr_start,
		unsigned long idxArr_end,
		unsigned translation_num,
		XFLOAT *weights,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		unsigned long current_oversampling)
{

	for (long unsigned i = 0; i < orientation_num*translation_num; i++)
		mapped_weights[i] = -999.;

	for (long unsigned i = idxArr_start; i < idxArr_end; i++)
		mapped_weights[ (rot_idx[i]-orientation_start) * translation_num + trans_idx[i] ]= weights[i];
}

void buildCorrImage(MlOptimiser *baseMLO, OptimisationParamters &op, CudaGlobalPtr<XFLOAT> &corr_img, long int ipart, long int group_id)
{
	// CC or not
	if((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			corr_img[n] = 1. / (op.local_sqrtXi2[ipart]*op.local_sqrtXi2[ipart]);
	else
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			corr_img[n] = *(op.local_Minvsigma2s[ipart].data + n );

	// ctf-correction or not ( NOTE this is not were the difference metric is ctf-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			corr_img[n] *= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
	// scale-correction or not ( NOTE this is not were the difference metric is scale-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
	if (baseMLO->do_scale_correction)
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			corr_img[n] *= myscale * myscale;
}

void generateEulerMatrices(
		XFLOAT padding_factor,
		ProjectionParams &ProjectionData,
		XFLOAT *eulers,
		bool inverse)
{
	double alpha, beta, gamma;
    double ca, sa, cb, sb, cg, sg;
    double cc, cs, sc, ss;

	for (long int i = 0; i < ProjectionData.rots.size(); i++)
	{
	    //TODO In a sense we're doing RAD2DEG just to do DEG2RAD here.
	    //The only place the degree value is actually used is in the metadata assignment.

	    alpha = DEG2RAD(ProjectionData.rots[i]);
	    beta  = DEG2RAD(ProjectionData.tilts[i]);
	    gamma = DEG2RAD(ProjectionData.psis[i]);

	    sincos(alpha, &sa, &ca);
	    sincos(beta,  &sb, &cb);
	    sincos(gamma, &sg, &cg);

	    cc = cb * ca;
	    cs = cb * sa;
	    sc = sb * ca;
	    ss = sb * sa;

		if(inverse)
		{
		    eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
		    eulers[9 * i + 1] = (-sg * cc - cg * sa) ;// * padding_factor; //10
		    eulers[9 * i + 2] = ( sc )               ;// * padding_factor; //20
		    eulers[9 * i + 3] = ( cg * cs + sg * ca) ;// * padding_factor; //01
		    eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
		    eulers[9 * i + 5] = ( ss )               ;// * padding_factor; //21
		    eulers[9 * i + 6] = (-cg * sb )          ;// * padding_factor; //02
		    eulers[9 * i + 7] = ( sg * sb )          ;// * padding_factor; //12
		    eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
		}
		else
		{
		    eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
		    eulers[9 * i + 1] = ( cg * cs + sg * ca) ;// * padding_factor; //01
		    eulers[9 * i + 2] = (-cg * sb )          ;// * padding_factor; //02
		    eulers[9 * i + 3] = (-sg * cc - cg * sa) ;// * padding_factor; //10
		    eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
		    eulers[9 * i + 5] = ( sg * sb )          ;// * padding_factor; //12
		    eulers[9 * i + 6] = ( sc )               ;// * padding_factor; //20
		    eulers[9 * i + 7] = ( ss )               ;// * padding_factor; //21
		    eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
		}
	}
}


long unsigned generateProjectionSetup(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		bool coarse,
		unsigned iclass,
		ProjectionParams &ProjectionData)
// FIXME : For coarse iteration this is **SLOW**    HERE ARE SOME NOTES FOR PARALLELIZING IT (GPU OFFLOAD):
/*
 *    Since it is based on push_back, parallelizing sould be fine given som atomic opreation appends,
 *    what takes time is looping through all this. The job-splitting in collect2jobs-preproccesing and
 *    divideOrientationsIntoBlockjobs() relies on chunks of shared orientations being adjacent in
 *    ProjectionData.rot_id (and thus also .rot_idx), but does not care which order those chunks appear
 *    in. So as long as a parallelilsm and "atomic push_back" is organised to use an orientation as a
 *    minimum unit, the job-splitting should be fine with the output.
 */
{
	//Local variables
	std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
	long int orientation_num = 0;

	for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
	{
		for (long int ipsi = sp.ipsi_min, ipart = 0; ipsi <= sp.ipsi_max; ipsi++, iorient++)
		{
			long int iorientclass = iclass * sp.nr_dir * sp.nr_psi + iorient;

			// Get prior for this direction and skip calculation if prior==0
			double pdf_orientation;
			if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
			{
				pdf_orientation = baseMLO->mymodel.pdf_class[iclass];
			}
			else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
			{
				pdf_orientation = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[iclass], idir);
			}
			else
			{
				pdf_orientation = op.directions_prior[idir] * op.psi_prior[ipsi];
			}
			// In the first pass, always proceed
			// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
			// if so, proceed with projecting the reference in that direction

			bool do_proceed = coarse ? true :
					baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant);

			if (do_proceed && pdf_orientation > 0.)
			{
				// Now get the oversampled (rot, tilt, psi) triplets
				// This will be only the original (rot,tilt,psi) triplet in the first pass (sp.current_oversampling==0)
				baseMLO->sampling.getOrientations(idir, ipsi, sp.current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

				// Loop over all oversampled orientations (only a single one in the first pass)
				for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++, ipart++)
				{
					ProjectionData.pushBackAll(	(long unsigned)iclass,
												oversampled_rot[iover_rot],
											    oversampled_tilt[iover_rot],
											    oversampled_psi[iover_rot],
											    iorientclass,
											    iover_rot 					);
					orientation_num ++;
				}
			}
		}
	}
	ProjectionData.orientation_num[iclass]=orientation_num;
	return orientation_num;
}

void runWavgKernel(
		CudaProjectorKernel &projector,
		XFLOAT *eulers,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *Fimgs_nomask_real,
		XFLOAT *Fimgs_nomask_imag,
		XFLOAT *sorted_weights,
		XFLOAT *ctfs,
		XFLOAT *Minvsigma2s,
		XFLOAT *wdiff2s_parts,
		XFLOAT *wdiff2s_AA,
		XFLOAT *wdiff2s_XA,
		XFLOAT *wavgs_real,
		XFLOAT *wavgs_imag,
		XFLOAT *Fweights,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		long int ipart,
		int group_id,
		int exp_iclass,
		XFLOAT part_scale,
		cudaStream_t stream)
{
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	unsigned orientation_chunks = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);
	CUDA_CPU_TIC("splitblocks");
	dim3 block_dim = splitCudaBlocks(orientation_chunks,false);
	CUDA_CPU_TOC("splitblocks");

	CUDA_CPU_TIC("cuda_kernel_wavg");

	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);

	if(projector.mdlZ!=0)
		cuda_kernel_wavg<true><<<block_dim,WAVG_BLOCK_SIZE,0,stream>>>(
			eulers,
			projector,
			image_size,
			orientation_num,
			Fimgs_real,
			Fimgs_imag,
			Fimgs_nomask_real,
			Fimgs_nomask_imag,
			sorted_weights,
			ctfs,
			Minvsigma2s,
			wdiff2s_parts,
			wdiff2s_AA,
			wdiff2s_XA,
			wavgs_real,
			wavgs_imag,
			Fweights,
			translation_num,
			(XFLOAT) op.sum_weight[ipart],
			(XFLOAT) op.significant_weight[ipart],
			baseMLO->refs_are_ctf_corrected,
			part_scale
			);
	else
		cuda_kernel_wavg<false><<<block_dim,WAVG_BLOCK_SIZE,0,stream>>>(
			eulers,
			projector,
			image_size,
			orientation_num,
			Fimgs_real,
			Fimgs_imag,
			Fimgs_nomask_real,
			Fimgs_nomask_imag,
			sorted_weights,
			ctfs,
			Minvsigma2s,
			wdiff2s_parts,
			wdiff2s_AA,
			wdiff2s_XA,
			wavgs_real,
			wavgs_imag,
			Fweights,
			translation_num,
			(XFLOAT) op.sum_weight[ipart],
			(XFLOAT) op.significant_weight[ipart],
			baseMLO->refs_are_ctf_corrected,
			part_scale
			);
	CUDA_CPU_TOC("cuda_kernel_wavg");
}

#define INIT_VALUE_BLOCK_SIZE 512
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

template< typename T>
void deviceInitValue(CudaGlobalPtr<T> data, T value)
{
	int grid_size = ceil((float)data.getSize()/(float)INIT_VALUE_BLOCK_SIZE);
	cuda_kernel_init_value<T><<< grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream() >>>(
			~data,
			value,
			data.getSize());
}



__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights
		)
{
	d_mweights[d_iorient[blockIdx.x] * blockDim.x + threadIdx.x] = d_allweights[blockIdx.x * blockDim.x + threadIdx.x];
}

void mapAllWeightsToMweights(
		unsigned long * d_iorient, //projectorPlan.iorientclasses
		XFLOAT * d_allweights, //allWeights
		XFLOAT * d_mweights, //Mweight
		unsigned long orientation_num, //projectorPlan.orientation_num
		unsigned long translation_num, //translation_num
		cudaStream_t stream
		)
{
	cuda_kernel_allweights_to_mweights<<< orientation_num, translation_num, 0, stream >>>(
			d_iorient,
			d_allweights,
			d_mweights);
}

void runDiff2KernelCoarse(
		CudaProjectorKernel &projector,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *d_eulers,
		XFLOAT *diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		int translation_num,
		int image_size,
		int ipart,
		int group_id,
		int exp_iclass,
		MlOptimiserCuda *cudaMLO,
		bool do_CC)
{
	if(!do_CC)
	{
		if(projector.mdlZ!=0)
		{
			if (translation_num > D2C_BLOCK_SIZE_3D)
			{
				printf("Number of coarse translations larger than %d not supported yet.\n", D2C_BLOCK_SIZE_3D);
				fflush(stdout);
				exit(1);
			}

			unsigned rest = orientation_num % D2C_EULERS_PER_BLOCK_3D;
			long unsigned even_orientation_num = orientation_num - rest;

			if (even_orientation_num != 0)
			{
				cuda_kernel_diff2_coarse<true, D2C_BLOCK_SIZE_3D, D2C_EULERS_PER_BLOCK_3D, 4>
				<<<even_orientation_num/D2C_EULERS_PER_BLOCK_3D,D2C_BLOCK_SIZE_3D,0,cudaMLO->classStreams[exp_iclass]>>>(
					d_eulers,
					trans_x,
					trans_y,
					Fimgs_real,
					Fimgs_imag,
					projector,
					corr_img,
					diff2s,
					translation_num,
					image_size);
			}

			if (rest != 0)
			{
				cuda_kernel_diff2_coarse<true, D2C_BLOCK_SIZE_3D, 1, 4>
				<<<rest,D2C_BLOCK_SIZE_3D,0,cudaMLO->classStreams[exp_iclass]>>>(
					&d_eulers[9*even_orientation_num],
					trans_x,
					trans_y,
					Fimgs_real,
					Fimgs_imag,
					projector,
					corr_img,
					&diff2s[translation_num*even_orientation_num],
					translation_num,
					image_size);
			}

		}
		else
		{

			if (translation_num > D2C_BLOCK_SIZE_2D)
			{
				printf("Number of coarse translations larger than %d not supported yet.\n", D2C_BLOCK_SIZE_2D);
				fflush(stdout);
				exit(1);
			}


			unsigned rest = orientation_num % D2C_EULERS_PER_BLOCK_2D;
			long unsigned even_orientation_num = orientation_num - rest;

			if (even_orientation_num != 0)
			{
				cuda_kernel_diff2_coarse<false, D2C_BLOCK_SIZE_2D, D2C_EULERS_PER_BLOCK_2D, 2>
				<<<even_orientation_num/D2C_EULERS_PER_BLOCK_2D,D2C_BLOCK_SIZE_2D,0,cudaMLO->classStreams[exp_iclass]>>>(
					d_eulers,
					trans_x,
					trans_y,
					Fimgs_real,
					Fimgs_imag,
					projector,
					corr_img,
					diff2s,
					translation_num,
					image_size);
			}

			if (rest != 0)
			{
				cuda_kernel_diff2_coarse<false, D2C_BLOCK_SIZE_2D, 1, 2>
				<<<rest,D2C_BLOCK_SIZE_2D,0,cudaMLO->classStreams[exp_iclass]>>>(
					&d_eulers[9*even_orientation_num],
					trans_x,
					trans_y,
					Fimgs_real,
					Fimgs_imag,
					projector,
					corr_img,
					&diff2s[translation_num*even_orientation_num],
					translation_num,
					image_size);
			}
		}
	}
	else
	{
		if(projector.mdlZ!=0)
			cuda_kernel_diff2_CC_coarse<true><<<orientation_num,BLOCK_SIZE,2*translation_num*BLOCK_SIZE*sizeof(XFLOAT),cudaMLO->classStreams[exp_iclass]>>>(
				d_eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,
				diff2s,
				translation_num,
				image_size,
				(XFLOAT) op.local_sqrtXi2[ipart]);
		else
			cuda_kernel_diff2_CC_coarse<false><<<orientation_num,BLOCK_SIZE,2*translation_num*BLOCK_SIZE*sizeof(XFLOAT),cudaMLO->classStreams[exp_iclass]>>>(
				d_eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,
				diff2s,
				translation_num,
				image_size,
				(XFLOAT) op.local_sqrtXi2[ipart]);
	}
}


void runDiff2KernelFine(
		CudaProjectorKernel &projector,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *eulers,
		long unsigned *rot_id,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		long unsigned *job_idx,
		long unsigned *job_num,
		XFLOAT *diff2s,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		long unsigned significant_num,
		unsigned image_size,
		int ipart,
		int exp_iclass,
		MlOptimiserCuda *cudaMLO,
		long unsigned job_num_count,
		bool do_CC)
{
    dim3 block_dim = splitCudaBlocks(job_num_count,false);

    if(!do_CC)
    {
		if(projector.mdlZ!=0)
			cuda_kernel_diff2_fine<true><<<block_dim,BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
				diff2s,
				image_size,
				op.highres_Xi2_imgs[ipart] / 2.,
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num);
		else
			cuda_kernel_diff2_fine<false><<<block_dim,BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
				diff2s,
				image_size,
				op.highres_Xi2_imgs[ipart] / 2.,
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num);
    }
    else
    {
		if(projector.mdlZ!=0)
			cuda_kernel_diff2_CC_fine<true><<<block_dim,BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,
				diff2s,
				image_size,
				op.highres_Xi2_imgs[ipart] / 2.,
				(XFLOAT) op.local_sqrtXi2[ipart],
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num);
		else
			cuda_kernel_diff2_CC_fine<false><<<block_dim,BLOCK_SIZE,0,cudaMLO->classStreams[exp_iclass]>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				projector,
				corr_img,
				diff2s,
				image_size,
				op.highres_Xi2_imgs[ipart] / 2.,
				(XFLOAT) op.local_sqrtXi2[ipart],
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num);
    }
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

	g_out_real[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp] = g_in_real[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp];
	g_out_imag[(kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp] = g_in_imag[(kp < 0 ? kp + iZ : kp)*iYX + (ip < 0 ? ip + iY : ip)*iX + jp];
}

void windowFourierTransform2(
		XFLOAT *d_in_real,
		XFLOAT *d_in_imag,
		XFLOAT *d_out_real,
		XFLOAT *d_out_imag,
		unsigned iX, unsigned iY, unsigned iZ, //Input dimensions
		unsigned oX, unsigned oY, unsigned oZ  //Output dimensions
		)
{
	if (iX > 1 && iY/2 + 1 != iX)
		REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");

	if (oY == iX)
		REPORT_ERROR("windowFourierTransform ERROR: there is a one-to-one map between input and output!");

	if (oY > iX)
	{
		long int max_r2 = (iX - 1) * (iX - 1);

		unsigned grid_dim = ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE);
		cuda_kernel_window_fourier_transform<true><<< grid_dim, WINDOW_FT_BLOCK_SIZE >>>(
				d_in_real,
				d_in_imag,
				d_out_real,
				d_out_imag,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				iX*iY*iZ,
				max_r2 );
	}
	else
	{
		unsigned grid_dim = ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE);
		cuda_kernel_window_fourier_transform<false><<< grid_dim, WINDOW_FT_BLOCK_SIZE >>>(
				d_in_real,
				d_in_imag,
				d_out_real,
				d_out_imag,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				oX*oY*oZ);
	}
}



void selfApplyBeamTilt2(MultidimArray<Complex > &Fimg, double beamtilt_x, double beamtilt_y,
		double wavelength, double Cs, double angpix, int ori_size)
{
	if (Fimg.getDim() != 2)
		REPORT_ERROR("applyBeamTilt can only be done on 2D Fourier Transforms!");

	double boxsize = angpix * ori_size;
	double factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize);

	for (unsigned n = 0 ; n < Fimg.yxdim; n ++)
	{
		unsigned i = n / Fimg.xdim;
		unsigned j = n % Fimg.xdim;
		unsigned jp = j;
		int ip = i < Fimg.xdim ? i : i - Fimg.ydim;

		double delta_phase = factor * (ip * ip + jp * jp) * (ip * beamtilt_y + jp * beamtilt_x);
		double realval = Fimg.data[i*Fimg.xdim+j].real;
		double imagval = Fimg.data[i*Fimg.xdim+j].imag;
		double mag = sqrt(realval * realval + imagval * imagval);
		double phas = atan2(imagval, realval) + DEG2RAD(delta_phase); // apply phase shift!
		realval = mag * cos(phas);
		imagval = mag * sin(phas);
		Fimg.data[i*Fimg.xdim+j] = Complex(realval, imagval);

	}

}
