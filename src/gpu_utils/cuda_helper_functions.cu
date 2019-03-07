#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_helper_functions.cuh"
#include "src/gpu_utils/cuda_kernels/BP.cuh"
#include "src/macros.h"
#include "src/error.h"

long int makeJobsForDiff2Fine(
		OptimisationParamters &op,  SamplingParameters &sp,
		long int orientation_num, long int translation_num,
		ProjectionParams &FineProjectionData,
		std::vector< long unsigned > &iover_transes,
		std::vector< long unsigned > &ihiddens,
		long int nr_over_orient, long int nr_over_trans, int ipart,
		IndexedDataArray &FPW, // FPW=FinePassWeights
		IndexedDataArrayMask &dataMask,
		int chunk)
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

				if(tk>=chunk)
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
		for(int i = 0; i < corr_img.getSize(); i++)
			corr_img[i] = 1. / (op.local_sqrtXi2[ipart]*op.local_sqrtXi2[ipart]);
	else
		for(int i = 0; i < corr_img.getSize(); i++)
			corr_img[i] = *(op.local_Minvsigma2s[ipart].data + i );

	// ctf-correction or not ( NOTE this is not were the difference metric is ctf-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
		for(int i = 0; i < corr_img.getSize(); i++)
			corr_img[i] *= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], i)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], i);
	// scale-correction or not ( NOTE this is not were the difference metric is scale-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
	if (baseMLO->do_scale_correction)
		for(int i = 0; i < corr_img.getSize(); i++)
			corr_img[i] *= myscale * myscale;
}

void generateEulerMatrices(
		XFLOAT padding_factor,
		ProjectionParams &ProjectionData,
		XFLOAT *eulers,
		bool inverse)
{
	RFLOAT alpha, beta, gamma;
    RFLOAT ca, sa, cb, sb, cg, sg;
    RFLOAT cc, cs, sc, ss;

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


long unsigned generateProjectionSetupFine(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
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
	std::vector< RFLOAT > oversampled_rot, oversampled_tilt, oversampled_psi;
	long int orientation_num = 0;

	for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
	{
		for (long int ipsi = sp.ipsi_min, ipart = 0; ipsi <= sp.ipsi_max; ipsi++, iorient++)
		{
			long int iorientclass = iclass * sp.nr_dir * sp.nr_psi + iorient;

			if (baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
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
		XFLOAT *Fimg_real,
		XFLOAT *Fimg_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT *sorted_weights,
		XFLOAT *ctfs,
		XFLOAT *wdiff2s_parts,
		XFLOAT *wdiff2s_AA,
		XFLOAT *wdiff2s_XA,
		OptimisationParamters &op,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		long int ipart,
		int group_id,
		int exp_iclass,
		XFLOAT part_scale,
		bool refs_are_ctf_corrected,
		bool data_is_3D,
		cudaStream_t stream)
{
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	dim3 block_dim = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);

	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);

	if (refs_are_ctf_corrected)
	{
		if(data_is_3D)
			cuda_kernel_wavg<true,true,true,WAVG_BLOCK_SIZE_DATA3D><<<block_dim,WAVG_BLOCK_SIZE_DATA3D,(3*WAVG_BLOCK_SIZE_DATA3D+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
		else if (projector.mdlZ!=0)
			cuda_kernel_wavg<true,true,false,WAVG_BLOCK_SIZE><<<block_dim,WAVG_BLOCK_SIZE,(3*WAVG_BLOCK_SIZE+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
		else
			cuda_kernel_wavg<true,false,false,WAVG_BLOCK_SIZE><<<block_dim,WAVG_BLOCK_SIZE,(3*WAVG_BLOCK_SIZE+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
	}
	else
	{
		if(data_is_3D)
			cuda_kernel_wavg<false,true,true,WAVG_BLOCK_SIZE_DATA3D><<<block_dim,WAVG_BLOCK_SIZE_DATA3D,(3*WAVG_BLOCK_SIZE_DATA3D+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
		else if (projector.mdlZ!=0)
			cuda_kernel_wavg<false,true,false,WAVG_BLOCK_SIZE><<<block_dim,WAVG_BLOCK_SIZE,(3*WAVG_BLOCK_SIZE+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
		else
			cuda_kernel_wavg<false,false,false,WAVG_BLOCK_SIZE><<<block_dim,WAVG_BLOCK_SIZE,(3*WAVG_BLOCK_SIZE+9)*sizeof(XFLOAT),stream>>>(
				eulers,
				projector,
				image_size,
				orientation_num,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				sorted_weights,
				ctfs,
				wdiff2s_parts,
				wdiff2s_AA,
				wdiff2s_XA,
				translation_num,
				(XFLOAT) op.sum_weight[ipart],
				(XFLOAT) op.significant_weight[ipart],
				part_scale
				);
	}
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
}

void runBackProjectKernel(
		CudaBackprojector &BP,
		CudaProjectorKernel &projector,
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
		bool data_is_3D,
		bool do_sgd,
		cudaStream_t optStream)
{

	if(BP.mdlZ==1)
	{
		cuda_kernel_backproject2D<<<imageCount,BP_2D_BLOCK_SIZE,0,optStream>>>(
				d_img_real,
				d_img_imag,
				trans_x,
				trans_y,
				d_weights,
				d_Minvsigma2s,
				d_ctfs,
				translation_num,
				significant_weight,
				weight_norm,
				d_eulers,
				BP.d_mdlReal,
				BP.d_mdlImag,
				BP.d_mdlWeight,
				BP.maxR,
				BP.maxR2,
				BP.padding_factor,
				imgX,
				imgY,
				imgX*imgY,
				BP.mdlX,
				BP.mdlInitY);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
	else
	{
		if(do_sgd)
		{
			if(data_is_3D)
				cuda_kernel_backprojectSGD<true><<<imageCount,BP_DATA3D_BLOCK_SIZE,0,optStream>>>(
					projector,
					d_img_real,
					d_img_imag,
					trans_x,
					trans_y,
					trans_z,
					d_weights,
					d_Minvsigma2s,
					d_ctfs,
					translation_num,
					significant_weight,
					weight_norm,
					d_eulers,
					BP.d_mdlReal,
					BP.d_mdlImag,
					BP.d_mdlWeight,
					BP.maxR,
					BP.maxR2,
					BP.padding_factor,
					imgX,
					imgY,
					imgZ,
					imgX*imgY*imgZ,
					BP.mdlX,
					BP.mdlY,
					BP.mdlInitY,
					BP.mdlInitZ);
			else
				cuda_kernel_backprojectSGD<false><<<imageCount,BP_REF3D_BLOCK_SIZE,0,optStream>>>(
					projector,
					d_img_real,
					d_img_imag,
					trans_x,
					trans_y,
					trans_z,
					d_weights,
					d_Minvsigma2s,
					d_ctfs,
					translation_num,
					significant_weight,
					weight_norm,
					d_eulers,
					BP.d_mdlReal,
					BP.d_mdlImag,
					BP.d_mdlWeight,
					BP.maxR,
					BP.maxR2,
					BP.padding_factor,
					imgX,
					imgY,
					imgZ,
					imgX*imgY*imgZ,
					BP.mdlX,
					BP.mdlY,
					BP.mdlInitY,
					BP.mdlInitZ);
		}
		else
		{
			if(data_is_3D)
				cuda_kernel_backproject3D<true><<<imageCount,BP_DATA3D_BLOCK_SIZE,0,optStream>>>(
					d_img_real,
					d_img_imag,
					trans_x,
					trans_y,
					trans_z,
					d_weights,
					d_Minvsigma2s,
					d_ctfs,
					translation_num,
					significant_weight,
					weight_norm,
					d_eulers,
					BP.d_mdlReal,
					BP.d_mdlImag,
					BP.d_mdlWeight,
					BP.maxR,
					BP.maxR2,
					BP.padding_factor,
					imgX,
					imgY,
					imgZ,
					imgX*imgY*imgZ,
					BP.mdlX,
					BP.mdlY,
					BP.mdlInitY,
					BP.mdlInitZ);
			else
				cuda_kernel_backproject3D<false><<<imageCount,BP_REF3D_BLOCK_SIZE,0,optStream>>>(
					d_img_real,
					d_img_imag,
					trans_x,
					trans_y,
					trans_z,
					d_weights,
					d_Minvsigma2s,
					d_ctfs,
					translation_num,
					significant_weight,
					weight_norm,
					d_eulers,
					BP.d_mdlReal,
					BP.d_mdlImag,
					BP.d_mdlWeight,
					BP.maxR,
					BP.maxR2,
					BP.padding_factor,
					imgX,
					imgY,
					imgZ,
					imgX*imgY*imgZ,
					BP.mdlX,
					BP.mdlY,
					BP.mdlInitY,
					BP.mdlInitZ);
		}
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
}



__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num
		)
{
	size_t idx = blockIdx.x * WEIGHT_MAP_BLOCK_SIZE + threadIdx.x;
	if (idx < orientation_num*translation_num)
		d_mweights[d_iorient[idx/translation_num] * translation_num + idx%translation_num] =
				d_allweights[idx/translation_num * translation_num + idx%translation_num];
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
	int grid_size = ceil((float)(orientation_num*translation_num)/(float)WEIGHT_MAP_BLOCK_SIZE);
	cuda_kernel_allweights_to_mweights<<< grid_size, WEIGHT_MAP_BLOCK_SIZE, 0, stream >>>(
			d_iorient,
			d_allweights,
			d_mweights,
			orientation_num,
			translation_num);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
}


size_t findThresholdIdxInCumulativeSum(CudaGlobalPtr<XFLOAT> &data, XFLOAT threshold)
{
	int grid_size = ceil((float)(data.getSize()-1)/(float)FIND_IN_CUMULATIVE_BLOCK_SIZE);
	if(grid_size==0)
	{
		return(0);
	}
	else
	{
		CudaGlobalPtr<size_t >  idx(1, data.getStream(), data.getAllocator());
		idx[0] = 0;
		idx.put_on_device();



		cuda_kernel_find_threshold_idx_in_cumulative<<< grid_size, FIND_IN_CUMULATIVE_BLOCK_SIZE, 0, data.getStream() >>>(
				~data,
				threshold,
				data.getSize()-1,
				~idx);
		idx.cp_to_host();
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(data.getStream()));

		return idx[0];
	}
}

void runDiff2KernelCoarse(
		CudaProjectorKernel &projector,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
		XFLOAT *corr_img,
		XFLOAT *Fimg_real,
		XFLOAT *Fimg_imag,
		XFLOAT *d_eulers,
		XFLOAT *diff2s,
		XFLOAT local_sqrtXi2,
		long unsigned orientation_num,
		int translation_num,
		int image_size,
		cudaStream_t stream,
		bool do_CC,
		bool data_is_3D)
{
	const int blocks3D = (data_is_3D? D2C_BLOCK_SIZE_DATA3D : D2C_BLOCK_SIZE_REF3D);

	if(!do_CC)
	{
		if(projector.mdlZ!=0)
		{

#ifdef CUDA_DOUBLE_PRECISION
			if (translation_num > blocks3D*4)
				CRITICAL(ERR_TRANSLIM);
#else
			if (translation_num > blocks3D*8)
				CRITICAL(ERR_TRANSLIM);
#endif

			unsigned rest = orientation_num % blocks3D;
			long unsigned even_orientation_num = orientation_num - rest;

			if (translation_num <= blocks3D)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D, D2C_EULERS_PER_BLOCK_DATA3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_DATA3D,D2C_BLOCK_SIZE_DATA3D,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D, D2C_EULERS_PER_BLOCK_REF3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_REF3D,D2C_BLOCK_SIZE_REF3D,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_DATA3D,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_REF3D,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
				}
			}
			else if (translation_num <= blocks3D*2)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*2, D2C_EULERS_PER_BLOCK_DATA3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_DATA3D,D2C_BLOCK_SIZE_DATA3D*2,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*2, D2C_EULERS_PER_BLOCK_REF3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_REF3D,D2C_BLOCK_SIZE_REF3D*2,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);

				}

				if (rest != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true, true, D2C_BLOCK_SIZE_DATA3D*2, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_DATA3D*2,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*2, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_REF3D*2,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
				}
			}
			else if (translation_num <= blocks3D*4)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*4, D2C_EULERS_PER_BLOCK_DATA3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_DATA3D,D2C_BLOCK_SIZE_DATA3D*4,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*4, D2C_EULERS_PER_BLOCK_REF3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_REF3D,D2C_BLOCK_SIZE_REF3D*4,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*4, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_DATA3D*4,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*4, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_REF3D*4,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
				}
			}
#ifndef CUDA_DOUBLE_PRECISION
			else
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*8, D2C_EULERS_PER_BLOCK_DATA3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_DATA3D,D2C_BLOCK_SIZE_DATA3D*8,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*8, D2C_EULERS_PER_BLOCK_REF3D, 4>
						<<<even_orientation_num/D2C_EULERS_PER_BLOCK_REF3D,D2C_BLOCK_SIZE_REF3D*8,0,stream>>>(
							d_eulers,
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							diff2s,
							translation_num,
							image_size);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						cuda_kernel_diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*8, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_DATA3D*8,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
					else
						cuda_kernel_diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*8, 1, 4>
						<<<rest,D2C_BLOCK_SIZE_REF3D*8,0,stream>>>(
							&d_eulers[9*even_orientation_num],
							trans_x,
							trans_y,
							trans_z,
							Fimg_real,
							Fimg_imag,
							projector,
							corr_img,
							&diff2s[translation_num*even_orientation_num],
							translation_num,
							image_size);
				}
			}
#endif
		}
		else
		{

			if (translation_num > D2C_BLOCK_SIZE_2D)
			{
				printf("Number of coarse translations larger than %d on the GPU not supported.\n", D2C_BLOCK_SIZE_2D);
				fflush(stdout);
				exit(1);
			}


			unsigned rest = orientation_num % D2C_EULERS_PER_BLOCK_2D;
			long unsigned even_orientation_num = orientation_num - rest;

			if (even_orientation_num != 0)
			{
				if(data_is_3D)
					cuda_kernel_diff2_coarse<false,true, D2C_BLOCK_SIZE_2D, D2C_EULERS_PER_BLOCK_2D, 2>
					<<<even_orientation_num/D2C_EULERS_PER_BLOCK_2D,D2C_BLOCK_SIZE_2D,0,stream>>>(
						d_eulers,
						trans_x,
						trans_y,
						trans_z,
						Fimg_real,
						Fimg_imag,
						projector,
						corr_img,
						diff2s,
						translation_num,
						image_size);
				else
					cuda_kernel_diff2_coarse<false,false, D2C_BLOCK_SIZE_2D, D2C_EULERS_PER_BLOCK_2D, 2>
					<<<even_orientation_num/D2C_EULERS_PER_BLOCK_2D,D2C_BLOCK_SIZE_2D,0,stream>>>(
						d_eulers,
						trans_x,
						trans_y,
						trans_z,
						Fimg_real,
						Fimg_imag,
						projector,
						corr_img,
						diff2s,
						translation_num,
						image_size);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
			}

			if (rest != 0)
			{
				if(data_is_3D)
					cuda_kernel_diff2_coarse<false,true, D2C_BLOCK_SIZE_2D, 1, 2>
					<<<rest,D2C_BLOCK_SIZE_2D,0,stream>>>(
						&d_eulers[9*even_orientation_num],
						trans_x,
						trans_y,
						trans_z,
						Fimg_real,
						Fimg_imag,
						projector,
						corr_img,
						&diff2s[translation_num*even_orientation_num],
						translation_num,
						image_size);
				else
					cuda_kernel_diff2_coarse<false,false, D2C_BLOCK_SIZE_2D, 1, 2>
					<<<rest,D2C_BLOCK_SIZE_2D,0,stream>>>(
						&d_eulers[9*even_orientation_num],
						trans_x,
						trans_y,
						trans_z,
						Fimg_real,
						Fimg_imag,
						projector,
						corr_img,
						&diff2s[translation_num*even_orientation_num],
						translation_num,
						image_size);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
			}
		}
	}
	else
	{
		dim3 CCblocks(orientation_num,translation_num);
		if(data_is_3D)
			cuda_kernel_diff2_CC_coarse<true,true,D2C_BLOCK_SIZE_DATA3D>
		<<<CCblocks,D2C_BLOCK_SIZE_DATA3D,0,stream>>>(
				d_eulers,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				projector,
				corr_img,
				diff2s,
				translation_num,
				image_size,
				local_sqrtXi2);
		else if(projector.mdlZ!=0)
			cuda_kernel_diff2_CC_coarse<true,false,D2C_BLOCK_SIZE_REF3D>
		<<<CCblocks,D2C_BLOCK_SIZE_REF3D,0,stream>>>(
				d_eulers,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				projector,
				corr_img,
				diff2s,
				translation_num,
				image_size,
				local_sqrtXi2);
		else
			cuda_kernel_diff2_CC_coarse<false,false,D2C_BLOCK_SIZE_2D>
		<<<CCblocks,D2C_BLOCK_SIZE_2D,0,stream>>>(
				d_eulers,
				Fimg_real,
				Fimg_imag,
				trans_x,
				trans_y,
				trans_z,
				projector,
				corr_img,
				diff2s,
				translation_num,
				image_size,
				local_sqrtXi2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
}


void runDiff2KernelFine(
		CudaProjectorKernel &projector,
		XFLOAT *corr_img,
		XFLOAT *Fimgs_real,
		XFLOAT *Fimgs_imag,
		XFLOAT *trans_x,
		XFLOAT *trans_y,
		XFLOAT *trans_z,
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
		cudaStream_t stream,
		long unsigned job_num_count,
		bool do_CC,
		bool data_is_3D)
{
    dim3 block_dim = job_num_count;

    if(!do_CC)
    {
			if(data_is_3D)
				cuda_kernel_diff2_fine<true,true, D2F_BLOCK_SIZE_DATA3D, D2F_CHUNK_DATA3D>
				<<<block_dim,D2F_BLOCK_SIZE_DATA3D,0,stream>>>(
					eulers,
					Fimgs_real,
					Fimgs_imag,
					trans_x,
					trans_y,
					trans_z,
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
			else if(projector.mdlZ!=0)
				cuda_kernel_diff2_fine<true,false,D2F_BLOCK_SIZE_REF3D,D2F_CHUNK_REF3D>
				<<<block_dim,D2F_BLOCK_SIZE_REF3D,0,stream>>>(
					eulers,
					Fimgs_real,
					Fimgs_imag,
					trans_x,
					trans_y,
					trans_z,
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
				cuda_kernel_diff2_fine<false,false,D2F_BLOCK_SIZE_2D,D2F_CHUNK_2D>
				<<<block_dim,D2F_BLOCK_SIZE_2D,0,stream>>>(
					eulers,
					Fimgs_real,
					Fimgs_imag,
					trans_x,
					trans_y,
					trans_z,
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
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }
    else
    {
    	if(data_is_3D)
			cuda_kernel_diff2_CC_fine<true,true,D2F_BLOCK_SIZE_DATA3D,D2F_CHUNK_DATA3D>
			<<<block_dim,D2F_BLOCK_SIZE_DATA3D,0,stream>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				trans_x,
				trans_y,
				trans_z,
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
    	else if(projector.mdlZ!=0)
			cuda_kernel_diff2_CC_fine<true,false,D2F_BLOCK_SIZE_REF3D,D2F_CHUNK_REF3D>
			<<<block_dim,D2F_BLOCK_SIZE_REF3D,0,stream>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				trans_x,
				trans_y,
				trans_z,
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
			cuda_kernel_diff2_CC_fine<false,false,D2F_BLOCK_SIZE_2D,D2F_CHUNK_2D>
			<<<block_dim,D2F_BLOCK_SIZE_2D,0,stream>>>(
				eulers,
				Fimgs_real,
				Fimgs_imag,
				trans_x,
				trans_y,
				trans_z,
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
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }

}

void runCollect2jobs(	dim3 grid_dim,
						XFLOAT * oo_otrans_x,          // otrans-size -> make const
						XFLOAT * oo_otrans_y,          // otrans-size -> make const
						XFLOAT * oo_otrans_z,          // otrans-size -> make const
						XFLOAT * myp_oo_otrans_x2y2z2, // otrans-size -> make const
						XFLOAT * weights,
						XFLOAT significant_weight,
						XFLOAT sum_weight,
						unsigned long nr_trans,
						unsigned long nr_oversampled_trans,
						unsigned long nr_oversampled_rot,
						int oversamples,
						bool skip_rots,
						XFLOAT * p_weights,
						XFLOAT * p_thr_wsum_prior_offsetx_class,
						XFLOAT * p_thr_wsum_prior_offsety_class,
						XFLOAT * p_thr_wsum_prior_offsetz_class,
						XFLOAT * p_thr_wsum_sigma2_offset,
						size_t * rot_idx,
						size_t * trans_idx,
						size_t * jobOrigin,
						size_t * jobExtent,
						bool data_is_3D
						)
{
	if(data_is_3D)
	{
		size_t shared_buffer = sizeof(XFLOAT)*SUMW_BLOCK_SIZE*5; // x+y+z+myp+weights
		cuda_kernel_collect2jobs<true><<<grid_dim,SUMW_BLOCK_SIZE,shared_buffer>>>(
				oo_otrans_x,          // otrans-size -> make const
				oo_otrans_y,          // otrans-size -> make const
				oo_otrans_z,          // otrans-size -> make const
				myp_oo_otrans_x2y2z2, // otrans-size -> make const
				weights,
				significant_weight,
				sum_weight,
				nr_trans,
				nr_oversampled_trans,
				nr_oversampled_rot,
				oversamples,
				skip_rots,
				p_weights,
				p_thr_wsum_prior_offsetx_class,
				p_thr_wsum_prior_offsety_class,
				p_thr_wsum_prior_offsetz_class,
				p_thr_wsum_sigma2_offset,
				rot_idx,
				trans_idx,
				jobOrigin,
				jobExtent);
	}
	else
	{
		size_t shared_buffer = sizeof(XFLOAT)*SUMW_BLOCK_SIZE*4; // x+y+myp+weights
		cuda_kernel_collect2jobs<false><<<grid_dim,SUMW_BLOCK_SIZE,shared_buffer>>>(
				oo_otrans_x,          // otrans-size -> make const
				oo_otrans_y,          // otrans-size -> make const
				oo_otrans_z,          // otrans-size -> make const
				myp_oo_otrans_x2y2z2, // otrans-size -> make const
				weights,
				significant_weight,
				sum_weight,
				nr_trans,
				nr_oversampled_trans,
				nr_oversampled_rot,
				oversamples,
				skip_rots,
				p_weights,
				p_thr_wsum_prior_offsetx_class,
				p_thr_wsum_prior_offsety_class,
				p_thr_wsum_prior_offsetz_class,
				p_thr_wsum_sigma2_offset,
				rot_idx,
				trans_idx,
				jobOrigin,
				jobExtent);
	}
}
//void windowFourierTransform2(
//		XFLOAT *d_in_real,
//		XFLOAT *d_in_imag,
//		XFLOAT *d_out_real,
//		XFLOAT *d_out_imag,
//		unsigned iX, unsigned iY, unsigned iZ, //Input dimensions
//		unsigned oX, unsigned oY, unsigned oZ,  //Output dimensions
//		cudaStream_t stream
//		)
//{
//	if (iX > 1 && iY/2 + 1 != iX)
//		REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
//
//	if (oY == iX)
//		REPORT_ERROR("windowFourierTransform ERROR: there is a one-to-one map between input and output!");
//
//	cudaMemInit<XFLOAT>( d_out_real, 0, (size_t) oX*oY*oZ, stream );
//	cudaMemInit<XFLOAT>( d_out_imag, 0, (size_t) oX*oY*oZ, stream );
//
//	if (oY > iX)
//	{
//		long int max_r2 = (iX - 1) * (iX - 1);
//
//		unsigned grid_dim = ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE);
//		cuda_kernel_window_fourier_transform<true><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, stream  >>>(
//				d_in_real,
//				d_in_imag,
//				d_out_real,
//				d_out_imag,
//				iX, iY, iZ, iX * iY, //Input dimensions
//				oX, oY, oZ, oX * oY, //Output dimensions
//				iX*iY*iZ,
//				max_r2 );
//	}
//	else
//	{
//		unsigned grid_dim = ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE);
//		cuda_kernel_window_fourier_transform<false><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, stream  >>>(
//				d_in_real,
//				d_in_imag,
//				d_out_real,
//				d_out_imag,
//				iX, iY, iZ, iX * iY, //Input dimensions
//				oX, oY, oZ, oX * oY, //Output dimensions
//				oX*oY*oZ);
//	}
//}

void windowFourierTransform2(
		CudaGlobalPtr<CUDACOMPLEX > &d_in,
		CudaGlobalPtr<CUDACOMPLEX > &d_out,
		size_t iX, size_t iY, size_t iZ, //Input dimensions
		size_t oX, size_t oY, size_t oZ,  //Output dimensions
		size_t Npsi,
		size_t pos,
		cudaStream_t stream)
{
	if (iX > 1 && iY/2 + 1 != iX)
		REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");

//	if (oX == iX)
//		REPORT_ERROR("windowFourierTransform ERROR: there is a one-to-one map between input and output!");


	deviceInitComplexValue(d_out, (XFLOAT)0.);
	HANDLE_ERROR(cudaStreamSynchronize(d_out.getStream()));

	if(oX==iX)
	{
		HANDLE_ERROR(cudaStreamSynchronize(d_in.getStream()));
		cudaCpyDeviceToDevice(&d_in.d_ptr[pos], ~d_out, oX*oY*oZ*Npsi, d_out.getStream() );
		return;
	}

	if (oX > iX)
	{
		long int max_r2 = (iX - 1) * (iX - 1);

		dim3 grid_dim(ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE),Npsi);
		cuda_kernel_window_fourier_transform<true><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, d_out.getStream() >>>(
				&d_in.d_ptr[pos],
				d_out.d_ptr,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				iX*iY*iZ,
				max_r2 );
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
	else
	{
		dim3 grid_dim(ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE),Npsi);
		cuda_kernel_window_fourier_transform<false><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, d_out.getStream() >>>(
				&d_in.d_ptr[pos],
				d_out.d_ptr,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				oX*oY*oZ);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
}



void selfApplyBeamTilt2(MultidimArray<Complex > &Fimg, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
		RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size)
{
	if (Fimg.getDim() != 2)
		REPORT_ERROR("applyBeamTilt can only be done on 2D Fourier Transforms!");

	RFLOAT boxsize = angpix * ori_size;
	RFLOAT factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize);

	for (unsigned n = 0 ; n < Fimg.yxdim; n ++)
	{
		unsigned i = n / Fimg.xdim;
		unsigned j = n % Fimg.xdim;
		unsigned jp = j;
		int ip = i < Fimg.xdim ? i : i - Fimg.ydim;

		RFLOAT delta_phase = factor * (ip * ip + jp * jp) * (ip * beamtilt_y + jp * beamtilt_x);
		RFLOAT realval = Fimg.data[i*Fimg.xdim+j].real;
		RFLOAT imagval = Fimg.data[i*Fimg.xdim+j].imag;
		RFLOAT mag = sqrt(realval * realval + imagval * imagval);
		RFLOAT phas = atan2(imagval, realval) + DEG2RAD(delta_phase); // apply phase shift!
		realval = mag * cos(phas);
		imagval = mag * sin(phas);
		Fimg.data[i*Fimg.xdim+j] = Complex(realval, imagval);

	}

}
