/*
#undef ALTCPU
#include <cuda_runtime.h>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_kernels/BP.cuh"
#include "src/macros.h"
#include "src/error.h"
*/

long int makeJobsForDiff2Fine(
		OptimisationParamters &op,  SamplingParameters &sp,
		long int orientation_num, long int translation_num,
		ProjectionParams &FineProjectionData,
		std::vector< long unsigned > &iover_transes,
		std::vector< long unsigned > &ihiddens,
		long int nr_over_orient, long int nr_over_trans, int img_id,
		IndexedDataArray &FPW, // FPW=FinePassWeights
		IndexedDataArrayMask &dataMask,
		int chunk)
{
	long unsigned w_base = dataMask.firstPos, w(0), k(0);
	// be on the safe side with the jobArrays: make them as large as they could possibly be
	// (this will be reduced at exit of this function)
	dataMask.setNumberOfJobs(orientation_num*translation_num);
	dataMask.setNumberOfWeights(orientation_num*translation_num);
	dataMask.jobOrigin.hostAlloc();
	dataMask.jobExtent.hostAlloc();

	dataMask.jobOrigin[k]=0;
	for (long unsigned i = 0; i < orientation_num; i++)
	{
		dataMask.jobExtent[k]=0;
		long int tk=0;
		long int iover_rot = FineProjectionData.iover_rots[i];
		for (long unsigned j = 0; j < translation_num; j++)
		{
			long int iover_trans = iover_transes[j];
			long int ihidden = FineProjectionData.iorientclasses[i] * sp.nr_trans + ihiddens[j];

			if(DIRECT_A2D_ELEM(op.Mcoarse_significant, img_id, ihidden)==1)
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

long int  makeJobsForCollect(IndexedDataArray &FPW,
    IndexedDataArrayMask &dataMask, unsigned long NewJobNum) // FPW=FinePassWeights
{
	// reset the old (diff2Fine) job-definitions
//	dataMask.jobOrigin.free_host();
//    dataMask.jobOrigin.free_device();
//    dataMask.jobExtent.free_host();
//    dataMask.jobExtent.free_device();
    dataMask.setNumberOfJobs(NewJobNum);
//    dataMask.jobOrigin.hostAlloc();
//    dataMask.jobExtent.hostAlloc();

	long int jobid=0;
	dataMask.jobOrigin[jobid]=0;
	dataMask.jobExtent[jobid]=1;
	long int crot =FPW.rot_idx[jobid]; // set current rot
	for(long int n=1; n<FPW.rot_idx.getSize(); n++)
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
//	dataMask.jobOrigin.putOnDevice();
//	dataMask.jobExtent.putOnDevice();

	return (jobid+1);
}

/*
 * Maps weights to a decoupled indexing of translations and orientations
 */
void mapWeights(
		unsigned long orientation_start,
		XFLOAT *mapped_weights,
		unsigned long orientation_num,
		unsigned long idxArr_start,
		unsigned long idxArr_end,
		unsigned long translation_num,
		XFLOAT *weights,
		long unsigned *rot_idx,
		long unsigned *trans_idx,
		unsigned long current_oversampling)
{

	for (long unsigned i = 0; i < orientation_num*translation_num; i++)
		mapped_weights[i] = -std::numeric_limits<XFLOAT>::max();

	for (long unsigned i = idxArr_start; i < idxArr_end; i++)
		mapped_weights[ (rot_idx[i]-orientation_start) * translation_num + trans_idx[i] ]= weights[i];
}

void buildCorrImage(MlOptimiser *baseMLO,
		OptimisationParamters &op,
		AccPtr<XFLOAT> &corr_img,
		int img_id,
		long int group_id)
{

	// CC or not
	if((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		for(size_t i = 0; i < corr_img.getSize(); i++)
			corr_img[i] = 1. / (op.local_sqrtXi2[img_id]*op.local_sqrtXi2[img_id]);
	else
		for(size_t i = 0; i < corr_img.getSize(); i++)
			corr_img[i] = *(op.local_Minvsigma2[img_id].data + i );

	// ctf-correction or not ( NOTE this is not were the difference metric is ctf-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	if (baseMLO->do_ctf_correction)
	{
		if (baseMLO->refs_are_ctf_corrected)
			for(size_t i = 0; i < corr_img.getSize(); i++)
				corr_img[i] *= DIRECT_MULTIDIM_ELEM(op.local_Fctf[img_id], i)*DIRECT_MULTIDIM_ELEM(op.local_Fctf[img_id], i);
	}

	// scale-correction or not ( NOTE this is not were the difference metric is scale-corrected, but
	// rather where we apply the additional correction to make the GPU-specific arithmetic equal
	// to the CPU method)
	XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
	if (baseMLO->do_scale_correction)
		for(size_t i = 0; i < corr_img.getSize(); i++)
			corr_img[i] *= myscale * myscale;
}

void generateEulerMatrices(
		ProjectionParams &ProjectionData,
		XFLOAT *eulers,
		bool inverse,
		Matrix2D<RFLOAT> &L,
		Matrix2D<RFLOAT> &R)
{
	RFLOAT alpha, beta, gamma;
    RFLOAT ca, sa, cb, sb, cg, sg;
    RFLOAT cc, cs, sc, ss;

    Matrix2D<RFLOAT> A(3,3);

    bool doL = (L.mdimx == 3 && L.mdimy == 3);
    bool doR = (R.mdimx == 3 && R.mdimy == 3);

	for (long int i = 0; i < ProjectionData.rots.size(); i++)
	{
	    //TODO In a sense we're doing RAD2DEG just to do DEG2RAD here.
	    //The only place the degree value is actually used is in the metadata assignment.

	    alpha = DEG2RAD(ProjectionData.rots[i]);
	    beta  = DEG2RAD(ProjectionData.tilts[i]);
	    gamma = DEG2RAD(ProjectionData.psis[i]);

#ifdef RELION_SINGLE_PRECISION
	    sincosf(alpha, &sa, &ca);
	    sincosf(beta,  &sb, &cb);
	    sincosf(gamma, &sg, &cg);
#else
	    sincos(alpha, &sa, &ca);
	    sincos(beta,  &sb, &cb);
	    sincos(gamma, &sg, &cg);
#endif

	    cc = cb * ca;
	    cs = cb * sa;
	    sc = sb * ca;
	    ss = sb * sa;

	    A(0, 0) =  cg * cc - sg * sa;
	    A(0, 1) =  cg * cs + sg * ca;
	    A(0, 2) = -cg * sb;
	    A(1, 0) = -sg * cc - cg * sa;
	    A(1, 1) = -sg * cs + cg * ca;
	    A(1, 2) = sg * sb;
	    A(2, 0) = sc;
	    A(2, 1) = ss;
	    A(2, 2) = cb;

		if (doL)
			A = L * A;

		if (doR)
			A = A * R;

		if(inverse)
			A = A.inv();

		for (int m = 0; m < 3; m ++)
			for (int n = 0; n < 3; n ++)
				eulers[9 * i + (m*3 + n)] = A(m, n);
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
		for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
		{
			long int iorientclass = iclass * sp.nr_dir * sp.nr_psi + iorient;

			if (baseMLO->isSignificantAnyImageAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
			{
				// Now get the oversampled (rot, tilt, psi) triplets
				// This will be only the original (rot,tilt,psi) triplet in the first pass (sp.current_oversampling==0)
				baseMLO->sampling.getOrientations(idir, ipsi, sp.current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

				// Loop over all oversampled orientations (only a single one in the first pass)
				for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
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
		AccProjectorKernel &projector,
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
		unsigned long image_size,
		int img_id,
		int group_id,
		int exp_iclass,
		XFLOAT part_scale,
		bool refs_are_ctf_corrected,
		bool ctf_premultiplied,
		bool data_is_3D,
		cudaStream_t stream)
{
	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);

	if (refs_are_ctf_corrected)
	{
		if(data_is_3D)
			AccUtilities::kernel_wavg<true,true,true,WAVG_BLOCK_SIZE_DATA3D>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
		else if (projector.mdlZ!=0)
			AccUtilities::kernel_wavg<true,true,false,WAVG_BLOCK_SIZE>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
		else
			AccUtilities::kernel_wavg<true,false,false,WAVG_BLOCK_SIZE>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
	}
	else
	{
		if(data_is_3D)
			AccUtilities::kernel_wavg<false,true,true,WAVG_BLOCK_SIZE_DATA3D>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
		else if (projector.mdlZ!=0)
			AccUtilities::kernel_wavg<false,true,false,WAVG_BLOCK_SIZE>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
		else
			AccUtilities::kernel_wavg<false,false,false,WAVG_BLOCK_SIZE>(
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
				(XFLOAT) op.sum_weight[img_id],
				(XFLOAT) op.significant_weight[img_id],
				part_scale,
				stream
				);
	}
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
}

void runBackProjectKernel(
		AccBackprojector &BP,
		AccProjectorKernel &projector,
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
		bool do_grad,
		bool ctf_premultiplied,
		cudaStream_t optStream)
{

	if(BP.mdlZ==1)
	{
#ifdef _CUDA_ENABLED
		if(do_grad)
			if(ctf_premultiplied)
				cuda_kernel_backproject2D_SGD<true><<<imageCount, BP_2D_BLOCK_SIZE, 0, optStream>>>(
						projector,
                                                d_img_real, d_img_imag,
						trans_x, trans_y,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgX*imgY,
						BP.mdlX, BP.mdlInitY);
			else
				cuda_kernel_backproject2D_SGD<false><<<imageCount, BP_2D_BLOCK_SIZE, 0, optStream>>>(
						projector,
                                                d_img_real, d_img_imag,
						trans_x, trans_y,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgX*imgY,
						BP.mdlX, BP.mdlInitY);
		else
			if(ctf_premultiplied)
				cuda_kernel_backproject2D<true><<<imageCount,BP_2D_BLOCK_SIZE,0,optStream>>>(
					d_img_real, d_img_imag,
					trans_x, trans_y,
					d_weights, d_Minvsigma2s, d_ctfs,
					translation_num, significant_weight, weight_norm, d_eulers,
					BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
					BP.maxR, BP.maxR2, BP.padding_factor,
					imgX, imgY, imgX*imgY,
					BP.mdlX, BP.mdlInitY);
			else
				cuda_kernel_backproject2D<false><<<imageCount,BP_2D_BLOCK_SIZE,0,optStream>>>(
					d_img_real, d_img_imag,
					trans_x, trans_y,
					d_weights, d_Minvsigma2s, d_ctfs,
					translation_num, significant_weight, weight_norm, d_eulers,
					BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
					BP.maxR, BP.maxR2, BP.padding_factor,
					imgX, imgY, imgX*imgY,
					BP.mdlX, BP.mdlInitY);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
		if(do_grad)
			if(ctf_premultiplied)
				CpuKernels::backproject2D_SGD<true>(
						imageCount, BP_2D_BLOCK_SIZE,
						projector,
                        d_img_real, d_img_imag,
                        trans_x, trans_y,
                        d_weights, d_Minvsigma2s, d_ctfs,
                        translation_num, significant_weight, weight_norm, d_eulers,
                        BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
                        BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
                        (unsigned)imgX, (unsigned)imgY, (unsigned)imgX*imgY,
                        (unsigned)BP.mdlX, BP.mdlInitY, BP.mutexes);
			else
				CpuKernels::backproject2D_SGD<false>(
						imageCount, BP_2D_BLOCK_SIZE,
						projector,
						d_img_real, d_img_imag,
						trans_x, trans_y,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgX*imgY,
						(unsigned)BP.mdlX, BP.mdlInitY, BP.mutexes);
		else
			if(ctf_premultiplied)
				CpuKernels::backproject2D<true>(
						imageCount, BP_2D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgX*imgY,
						(unsigned)BP.mdlX, BP.mdlInitY, BP.mutexes);
			else
				CpuKernels::backproject2D<false>(
						imageCount, BP_2D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgX*imgY,
						(unsigned)BP.mdlX, BP.mdlInitY, BP.mutexes);
#endif
	}
	else
	{
		if(do_grad)
		{
			if(data_is_3D)
#ifdef _CUDA_ENABLED
				if(ctf_premultiplied)
                    cuda_kernel_backproject3D_SGD<true, true><<<imageCount, BP_DATA3D_BLOCK_SIZE, 0, optStream>>>(
                            projector, d_img_real, d_img_imag,
                            trans_x, trans_y, trans_z,
                            d_weights, d_Minvsigma2s, d_ctfs,
                            translation_num, significant_weight, weight_norm, d_eulers,
                            BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
                            BP.maxR, BP.maxR2, BP.padding_factor,
                            imgX, imgY, imgZ, imgX * imgY * imgZ,
                            BP.mdlX, BP.mdlY, BP.mdlInitY, BP.mdlInitZ);
				else
                    cuda_kernel_backproject3D_SGD<true, false><<<imageCount, BP_DATA3D_BLOCK_SIZE, 0, optStream>>>(
                            projector, d_img_real, d_img_imag,
                            trans_x, trans_y, trans_z,
                            d_weights, d_Minvsigma2s, d_ctfs,
                            translation_num, significant_weight, weight_norm, d_eulers,
                            BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
                            BP.maxR, BP.maxR2, BP.padding_factor,
                            imgX, imgY, imgZ, imgX * imgY * imgZ,
                            BP.mdlX, BP.mdlY, BP.mdlInitY, BP.mdlInitZ);
#else
				if(ctf_premultiplied)
					CpuKernels::backproject3D_SGD<true, true>(imageCount, BP_DATA3D_BLOCK_SIZE,
						projector, d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
				else
					CpuKernels::backproject3D_SGD<true, false>(imageCount, BP_DATA3D_BLOCK_SIZE,
						projector, d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);

#endif
			else
#ifdef _CUDA_ENABLED
		        if(ctf_premultiplied)
				cuda_kernel_backproject3D_SGD<false, true><<<imageCount, BP_REF3D_BLOCK_SIZE, 0, optStream>>>(
			                projector, d_img_real, d_img_imag,
			                trans_x, trans_y, trans_z,
			                d_weights, d_Minvsigma2s, d_ctfs,
			                translation_num, significant_weight, weight_norm, d_eulers,
			                BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
			                BP.maxR, BP.maxR2, BP.padding_factor,
			                imgX, imgY, imgZ, imgX * imgY * imgZ,
			                BP.mdlX, BP.mdlY, BP.mdlInitY, BP.mdlInitZ);
			    else
			        cuda_kernel_backproject3D_SGD<false, false><<<imageCount, BP_REF3D_BLOCK_SIZE, 0, optStream>>>(
			                projector, d_img_real, d_img_imag,
			                trans_x, trans_y, trans_z,
			                d_weights, d_Minvsigma2s, d_ctfs,
			                translation_num, significant_weight, weight_norm, d_eulers,
			                BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
			                BP.maxR, BP.maxR2, BP.padding_factor,
			                imgX, imgY, imgZ, imgX * imgY * imgZ,
			                BP.mdlX, BP.mdlY, BP.mdlInitY, BP.mdlInitZ);

#else
				if(ctf_premultiplied)
					CpuKernels::backproject3D_SGD<false, true>(imageCount, BP_REF3D_BLOCK_SIZE,
						projector, d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
				else
					CpuKernels::backproject3D_SGD<false, false>(imageCount, BP_REF3D_BLOCK_SIZE,
						projector, d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
#endif
		}
		else
		{
			if(data_is_3D)
#ifdef _CUDA_ENABLED
				if(ctf_premultiplied)
					cuda_kernel_backproject3D<true, true><<<imageCount,BP_DATA3D_BLOCK_SIZE,0,optStream>>>(
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, imgX*imgY*imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ);
				else
					cuda_kernel_backproject3D<true, false><<<imageCount,BP_DATA3D_BLOCK_SIZE,0,optStream>>>(
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, imgX*imgY*imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ);

#else
			    if(ctf_premultiplied)
					CpuKernels::backproject3D<true, true>(imageCount,BP_DATA3D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
			    else
					CpuKernels::backproject3D<true, false>(imageCount,BP_DATA3D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);

#endif
			else
#ifdef _CUDA_ENABLED
			    if(ctf_premultiplied)
					cuda_kernel_backproject3D<false, true><<<imageCount,BP_REF3D_BLOCK_SIZE,0,optStream>>>(
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, imgX*imgY*imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ);
			    else
					cuda_kernel_backproject3D<false, false><<<imageCount,BP_REF3D_BLOCK_SIZE,0,optStream>>>(
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, BP.padding_factor,
						imgX, imgY, imgZ, imgX*imgY*imgZ,
						BP.mdlX, BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ);

#else
#if 1 //TODO Clean this up
			if(ctf_premultiplied)
				CpuKernels::backprojectRef3D<true>(imageCount,
					d_img_real, d_img_imag,
					trans_x, trans_y,
					d_weights, d_Minvsigma2s, d_ctfs,
					translation_num, significant_weight, weight_norm, d_eulers,
					BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
					BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
					(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
					(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
			else
				CpuKernels::backprojectRef3D<false>(imageCount,
					d_img_real, d_img_imag,
					trans_x, trans_y,
					d_weights, d_Minvsigma2s, d_ctfs,
					translation_num, significant_weight, weight_norm, d_eulers,
					BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
					BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
					(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
					(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);

#else
				if(ctf_premultiplied)
					CpuKernels::backproject3D<false, true>(imageCount,BP_REF3D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
				else
					CpuKernels::backproject3D<false, false>(imageCount,BP_REF3D_BLOCK_SIZE,
						d_img_real, d_img_imag,
						trans_x, trans_y, trans_z,
						d_weights, d_Minvsigma2s, d_ctfs,
						translation_num, significant_weight, weight_norm, d_eulers,
						BP.d_mdlReal, BP.d_mdlImag, BP.d_mdlWeight,
						BP.maxR, BP.maxR2, (XFLOAT)BP.padding_factor,
						(unsigned)imgX, (unsigned)imgY, (unsigned)imgZ, (size_t)imgX*(size_t)imgY*(size_t)imgZ,
						(unsigned)BP.mdlX, (unsigned)BP.mdlY, BP.mdlInitY, 	BP.mdlInitZ, BP.mutexes);
#endif
#endif
		} // do_grad is false
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	}
}

#define WEIGHT_MAP_BLOCK_SIZE 512

void mapAllWeightsToMweights(
		unsigned long * d_iorient, //projectorPlan.iorientclasses
		XFLOAT * d_allweights, //allWeights
		XFLOAT * d_mweights, //Mweight
		unsigned long orientation_num, //projectorPlan.orientation_num
		unsigned long translation_num, //translation_num
		cudaStream_t stream
		)
{
	size_t combinations = orientation_num*translation_num;
	int grid_size = ceil((float)(combinations)/(float)WEIGHT_MAP_BLOCK_SIZE);
#ifdef _CUDA_ENABLED
	cuda_kernel_allweights_to_mweights<<< grid_size, WEIGHT_MAP_BLOCK_SIZE, 0, stream >>>(
			d_iorient,
			d_allweights,
			d_mweights,
			orientation_num,
			translation_num,
			WEIGHT_MAP_BLOCK_SIZE);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
	for (size_t i=0; i < combinations; i++)
		d_mweights[d_iorient[i/translation_num] * translation_num + i%translation_num] =
			d_allweights[i/translation_num * translation_num + i%translation_num];
			// TODO - isn't this just d_allweights[idx + idx%translation_num]?   Really?
#endif
}

void runDiff2KernelCoarse(
		AccProjectorKernel &projector,
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
		long unsigned translation_num,
		long unsigned image_size,
		cudaStream_t stream,
		bool do_CC,
		bool data_is_3D)
{
	const long unsigned blocks3D = (data_is_3D? D2C_BLOCK_SIZE_DATA3D : D2C_BLOCK_SIZE_REF3D);

	if(!do_CC)
	{
		if(projector.mdlZ!=0)
		{

#ifdef ACC_DOUBLE_PRECISION
			if (translation_num > blocks3D*4)
				CRITICAL(ERR_TRANSLIM);
#else
			if (translation_num > blocks3D*8)
				CRITICAL(ERR_TRANSLIM);
#endif

			long unsigned rest = orientation_num % blocks3D;
			long unsigned even_orientation_num = orientation_num - rest;
// TODO - find a more compact way to represent these combinations resulting in
// a single call to diff2_course?
			if (translation_num <= blocks3D)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D, D2C_EULERS_PER_BLOCK_DATA3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_DATA3D,
							D2C_BLOCK_SIZE_DATA3D,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D, D2C_EULERS_PER_BLOCK_REF3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_REF3D,
							D2C_BLOCK_SIZE_REF3D,
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
							image_size,
							stream);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_DATA3D,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_REF3D,
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
							image_size,
							stream);
				}
			}
			else if (translation_num <= blocks3D*2)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*2, D2C_EULERS_PER_BLOCK_DATA3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_DATA3D,
							D2C_BLOCK_SIZE_DATA3D*2,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*2, D2C_EULERS_PER_BLOCK_REF3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_REF3D,
							D2C_BLOCK_SIZE_REF3D*2,
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
							image_size,
							stream);

				}

				if (rest != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true, true, D2C_BLOCK_SIZE_DATA3D*2, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_DATA3D*2,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*2, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_REF3D*2,
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
							image_size,
							stream);
				}
			}
			else if (translation_num <= blocks3D*4)
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*4, D2C_EULERS_PER_BLOCK_DATA3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_DATA3D,
							D2C_BLOCK_SIZE_DATA3D*4,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*4, D2C_EULERS_PER_BLOCK_REF3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_REF3D,
							D2C_BLOCK_SIZE_REF3D*4,
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
							image_size,
							stream);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*4, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_DATA3D*4,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*4, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_REF3D*4,
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
							image_size,
							stream);
				}
			}
#ifndef ACC_DOUBLE_PRECISION
			else
			{
				if (even_orientation_num != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*8, D2C_EULERS_PER_BLOCK_DATA3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_DATA3D,
							D2C_BLOCK_SIZE_DATA3D*8,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*8, D2C_EULERS_PER_BLOCK_REF3D, 4>(
							even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_REF3D,
							D2C_BLOCK_SIZE_REF3D*8,
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
							image_size,
							stream);
				}

				if (rest != 0)
				{
					if(data_is_3D)
						AccUtilities::diff2_coarse<true,true, D2C_BLOCK_SIZE_DATA3D*8, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_DATA3D*8,
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
							image_size,
							stream);
					else
						AccUtilities::diff2_coarse<true,false, D2C_BLOCK_SIZE_REF3D*8, 1, 4>(
							rest,
							D2C_BLOCK_SIZE_REF3D*8,
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
							image_size,
							stream);
				}
			}
#endif
		}  // projector.mdlZ!=0
		else
		{

			if (translation_num > D2C_BLOCK_SIZE_2D)
			{
				printf("Number of coarse translations larger than %d on the GPU not supported.\n", D2C_BLOCK_SIZE_2D);
				fflush(stdout);
				exit(1);
			}


			long unsigned rest = orientation_num % (unsigned long)D2C_EULERS_PER_BLOCK_2D;
			long unsigned even_orientation_num = orientation_num - rest;

			if (even_orientation_num != 0)
			{
				if(data_is_3D)
					AccUtilities::diff2_coarse<false,true, D2C_BLOCK_SIZE_2D, D2C_EULERS_PER_BLOCK_2D, 2>(
						even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_2D,
						D2C_BLOCK_SIZE_2D,
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
						image_size,
						stream);
				else
					AccUtilities::diff2_coarse<false,false, D2C_BLOCK_SIZE_2D, D2C_EULERS_PER_BLOCK_2D, 2>(
						even_orientation_num/(unsigned long)D2C_EULERS_PER_BLOCK_2D,
						D2C_BLOCK_SIZE_2D,
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
						image_size,
						stream);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
			}

			if (rest != 0)
			{
				if(data_is_3D)
					AccUtilities::diff2_coarse<false,true, D2C_BLOCK_SIZE_2D, 1, 2>(
						rest,
						D2C_BLOCK_SIZE_2D,
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
						image_size,
						stream);
				else
					AccUtilities::diff2_coarse<false,false, D2C_BLOCK_SIZE_2D, 1, 2>(
						rest,
						D2C_BLOCK_SIZE_2D,
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
						image_size,
						stream);
				LAUNCH_HANDLE_ERROR(cudaGetLastError());
			}
		}  // projector.mdlZ==0
	}  // !do_CC
	else
	{  // do_CC
// TODO - find a more compact way to represent these combinations resulting in
// a single call to diff2_CC_course?
		// dim3 CCblocks(orientation_num,translation_num);
		if(data_is_3D)
			AccUtilities::diff2_CC_coarse<true,true,D2C_BLOCK_SIZE_DATA3D>(
				orientation_num,
				D2C_BLOCK_SIZE_DATA3D,
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
				local_sqrtXi2,
				stream);
		else if(projector.mdlZ!=0)
			AccUtilities::diff2_CC_coarse<true,false,D2C_BLOCK_SIZE_REF3D>(
				orientation_num,
				D2C_BLOCK_SIZE_REF3D,
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
				local_sqrtXi2,
				stream);
		else
			AccUtilities::diff2_CC_coarse<false,false,D2C_BLOCK_SIZE_2D>(
				orientation_num,
				D2C_BLOCK_SIZE_2D,
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
				local_sqrtXi2,
				stream);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
	} // do_CC
}


void runDiff2KernelFine(
		AccProjectorKernel &projector,
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
		unsigned long image_size,
		int img_id,
		int exp_iclass,
		cudaStream_t stream,
		long unsigned job_num_count,
		bool do_CC,
		bool data_is_3D)
{
    long unsigned block_dim = job_num_count;

    if(!do_CC)
    {
// TODO - find a more compact way to represent these combinations resulting in
// a single call to diff2_fine?
			if(data_is_3D)
				AccUtilities::diff2_fine<true,true, D2F_BLOCK_SIZE_DATA3D, D2F_CHUNK_DATA3D>(
					block_dim,
					D2F_BLOCK_SIZE_DATA3D,
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
					op.highres_Xi2_img[img_id] / 2.,
					orientation_num,
					translation_num,
					job_num_count, //significant_num,
					rot_idx,
					trans_idx,
					job_idx,
					job_num,
					stream);
			else if(projector.mdlZ!=0)
				AccUtilities::diff2_fine<true,false,D2F_BLOCK_SIZE_REF3D,D2F_CHUNK_REF3D>(
					block_dim,
					D2F_BLOCK_SIZE_REF3D,
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
					op.highres_Xi2_img[img_id] / 2.,
					orientation_num,
					translation_num,
					job_num_count, //significant_num,
					rot_idx,
					trans_idx,
					job_idx,
					job_num,
					stream);
			else
				AccUtilities::diff2_fine<false,false,D2F_BLOCK_SIZE_2D,D2F_CHUNK_2D>(
					block_dim,
					D2F_BLOCK_SIZE_2D,
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
					op.highres_Xi2_img[img_id] / 2.,
					orientation_num,
					translation_num,
					job_num_count, //significant_num,
					rot_idx,
					trans_idx,
					job_idx,
					job_num,
					stream);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }
    else
    {
// TODO - find a more compact way to represent these combinations resulting in
// a single call to diff2_CC_fine?
    	if(data_is_3D)
			AccUtilities::diff2_CC_fine<true,true,D2F_BLOCK_SIZE_DATA3D,D2F_CHUNK_DATA3D>(
				block_dim,
				D2F_BLOCK_SIZE_DATA3D,
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
				op.highres_Xi2_img[img_id] / 2.,
				(XFLOAT) op.local_sqrtXi2[img_id],
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num,
				stream);
    	else if(projector.mdlZ!=0)
			AccUtilities::diff2_CC_fine<true,false,D2F_BLOCK_SIZE_REF3D,D2F_CHUNK_REF3D>(
				block_dim,
				D2F_BLOCK_SIZE_REF3D,
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
				op.highres_Xi2_img[img_id] / 2.,
				(XFLOAT) op.local_sqrtXi2[img_id],
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num,
				stream);
		else
			AccUtilities::diff2_CC_fine<false,false,D2F_BLOCK_SIZE_2D,D2F_CHUNK_2D>(
				block_dim,
				D2F_BLOCK_SIZE_2D,
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
				op.highres_Xi2_img[img_id] / 2.,
				(XFLOAT) op.local_sqrtXi2[img_id],
				orientation_num,
				translation_num,
				job_num_count, //significant_num,
				rot_idx,
				trans_idx,
				job_idx,
				job_num,
				stream);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }

}

void runCollect2jobs(	int grid_dim,
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
						unsigned long oversamples,
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
	if (data_is_3D) {
#ifdef _CUDA_ENABLED
	dim3 numblocks(grid_dim);
	size_t shared_buffer = sizeof(XFLOAT)*SUMW_BLOCK_SIZE*5; // x+y+z+myp+weights
	cuda_kernel_collect2jobs<true><<<numblocks,SUMW_BLOCK_SIZE,shared_buffer>>>(
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
#else
		CpuKernels::collect2jobs<true>(grid_dim, SUMW_BLOCK_SIZE,
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
#endif
	}
	else
	{
#ifdef _CUDA_ENABLED
	dim3 numblocks(grid_dim);
	size_t shared_buffer = sizeof(XFLOAT)*SUMW_BLOCK_SIZE*4; // x+y+myp+weights
	cuda_kernel_collect2jobs<false><<<numblocks,SUMW_BLOCK_SIZE,shared_buffer>>>(
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
#else
	CpuKernels::collect2jobs<false>(grid_dim, SUMW_BLOCK_SIZE,
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
#endif
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
#define WINDOW_FT_BLOCK_SIZE 128

void windowFourierTransform2(
		AccPtr<ACCCOMPLEX > &d_in,
		AccPtr<ACCCOMPLEX > &d_out,
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


	deviceInitComplexValue<ACCCOMPLEX>(d_out, (XFLOAT)0.);
	HANDLE_ERROR(cudaStreamSynchronize(d_out.getStream()));

	if(oX==iX)
	{
		HANDLE_ERROR(cudaStreamSynchronize(d_in.getStream()));
#ifdef _CUDA_ENABLED
		cudaCpyDeviceToDevice(&d_in(pos), ~d_out, oX*oY*oZ*Npsi, d_out.getStream() );
#else
		memcpy(&d_out[0], &d_in[0], oX*oY*oZ*Npsi*sizeof(ACCCOMPLEX));
#endif
		return;
	}

	if (oX > iX)
	{
		long int max_r2 = (iX - 1) * (iX - 1);

#ifdef _CUDA_ENABLED
		dim3 grid_dim(ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE),Npsi);
		cuda_kernel_window_fourier_transform<true><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, d_out.getStream() >>>(
				&d_in(pos),
				~d_out,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				iX*iY*iZ,
				WINDOW_FT_BLOCK_SIZE,
				max_r2);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
		size_t grid_dim = (size_t)( ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE));
		CpuKernels::window_fourier_transform<true>(
				grid_dim,
				Npsi,
				WINDOW_FT_BLOCK_SIZE,
				&d_in[pos],
				&d_out[0],
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				iX*iY*iZ,
				max_r2);
#endif
	}
	else
	{
#ifdef _CUDA_ENABLED
		dim3 grid_dim(ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE),Npsi);
		cuda_kernel_window_fourier_transform<false><<< grid_dim, WINDOW_FT_BLOCK_SIZE, 0, d_out.getStream() >>>(
				&d_in(pos),
				~d_out,
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				oX*oY*oZ,
				WINDOW_FT_BLOCK_SIZE);
		LAUNCH_HANDLE_ERROR(cudaGetLastError());
#else
		int grid_dim = (int)( ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE));
		CpuKernels::window_fourier_transform<false>(
				grid_dim,
				Npsi,
				WINDOW_FT_BLOCK_SIZE,
				&d_in[pos],
				&d_out[0],
				iX, iY, iZ, iX * iY, //Input dimensions
				oX, oY, oZ, oX * oY, //Output dimensions
				oX*oY*oZ
				);
#endif
	}
}

void run_calcPowerSpectrum(Complex *dFaux, int padoridim, Complex *ddata, int data_sz, RFLOAT *dpower_spectrum, RFLOAT *dcounter,
											  int max_r2, int min_r2, RFLOAT normfft, RFLOAT padding_factor, RFLOAT weight,
											  RFLOAT *dfourier_mask, int fx, int fy, int fz, bool do_fourier_mask, bool if3D)
{
#ifdef CUDA
	dim3 bs(32,4);
	dim3 gs(ceil((padoridim/2+1)/(float)bs.x), ceil(padoridim/(float)bs.y));
	if(if3D)
	{
		bs.z = 2;
		gs.z = ceil(padoridim/(float)bs.z); 
	}
	if(sizeof(RFLOAT) == sizeof(double))
		cuda_kernel_calcPowerSpectrum<<<gs,bs>>>((double2*)dFaux,padoridim,(double2*)ddata,data_sz,dpower_spectrum,dcounter,
												  max_r2,min_r2,normfft,padding_factor,weight,dfourier_mask,fx,fy,fz,do_fourier_mask);
	else
		cuda_kernel_calcPowerSpectrum<<<gs,bs>>>((float2*)dFaux,padoridim,(float2*)ddata,data_sz,dpower_spectrum,dcounter,
												  max_r2,min_r2,normfft,padding_factor,weight,dfourier_mask,fx,fy,fz,do_fourier_mask);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#endif
}

void run_updatePowerSpectrum(RFLOAT *dcounter, int sz, RFLOAT *dpower_spectrum)
{
#ifdef CUDA
	cuda_kernel_updatePowerSpectrum<<<ceil(sz/(float)256),256>>>(dcounter, dpower_spectrum, sz);
	LAUNCH_HANDLE_ERROR(cudaGetLastError());
#endif
}

void scale(RFLOAT *img, size_t sz, RFLOAT val, cudaStream_t stream)
{
	int block_size = 256;
	int MultiBsize = ceil(sz/(float)block_size);
#ifdef CUDA
	AccUtilities::multiply(MultiBsize,block_size, stream, img, val, (size_t)sz);
#endif
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
