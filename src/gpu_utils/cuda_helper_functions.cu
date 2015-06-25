#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <signal.h>


long int divideOrientationsIntoBlockjobs( OptimisationParamters &op,  SamplingParameters &sp,
										  long int orientation_num, long int translation_num,
										  std::vector< long unsigned > &iorientclasses,
										  std::vector< long unsigned > &iover_rots,
										  std::vector< long unsigned > &iover_transes,
										  std::vector< long unsigned > &ihiddens,
										  long int nr_over_orient, long int nr_over_trans, int ipart,
										  CudaGlobalPtr <long unsigned> &rot_id,
										  CudaGlobalPtr <long unsigned> &rot_idx,
										  CudaGlobalPtr <long unsigned> &trans_idx,
										  CudaGlobalPtr <long unsigned> &ihidden_overs,
										  CudaGlobalPtr <long unsigned> &job_idx,
										  CudaGlobalPtr <long unsigned> &job_num)
{
	rot_id.size=		orientation_num*translation_num;
	rot_idx.size=		orientation_num*translation_num;
	trans_idx.size=		orientation_num*translation_num;
	ihidden_overs.size=	orientation_num*translation_num;
	job_idx.size=		orientation_num*translation_num;
	job_num.size=		orientation_num*translation_num;

	rot_id.host_alloc();
	rot_idx.host_alloc();
	trans_idx.host_alloc();
	ihidden_overs.host_alloc();
	job_idx.host_alloc();
	job_num.host_alloc();


	long int significant_num(0), k(0);

	job_idx[k]=0;
	for (long unsigned i = 0; i < orientation_num; i++)
	{
		job_num[k]=0;
		int tk=0;
		long int iover_rot = iover_rots[i];
		for (long unsigned j = 0; j < translation_num; j++)
		{
			long int iover_trans = iover_transes[j];
			long int ihidden = iorientclasses[i] * sp.nr_trans + ihiddens[j];

			if(DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden)==1)
			{
				rot_id[significant_num] = iorientclasses[i]; 	// where to look for priors etc
				rot_idx[significant_num] = i;					// which rot for this significant task
				trans_idx[significant_num] = j;					// which trans       - || -
				ihidden_overs[significant_num]= (ihidden * nr_over_orient + iover_rot) * nr_over_trans + iover_trans;

				if(tk>=PROJDIFF_CHUNK_SIZE)
				{
					tk=0;             // reset counter
					k++;              // use new element
					job_idx[k]=significant_num;
					job_num[k]=0;   // prepare next element for ++ incrementing
				}
				tk++;                 // increment limit
				job_num[k]++;       // increment number of transes this ProjDiff-block
				significant_num++;
			}
			else if(tk!=0) // start a new one with the same rotidx - we expect transes to be sequential.
			{
				tk=0;             // reset counter
				k++;              // use new element
				job_idx[k]=significant_num;
				job_num[k]=0;   // prepare next element for ++ incrementing
			}
		}
		if(tk>0) // use new element (if tk==0) then we are currently on an element with no signif, so we should continue using this element
		{
			k++;
			job_idx[k]=significant_num;
			job_num[k]=0;
		}
	}

	job_num.size		=k;
	job_idx.size		=k;
	rot_id.size 		=significant_num;
	rot_idx.size		=significant_num;
	trans_idx.size		=significant_num;
	ihidden_overs.size  =significant_num;

	return(significant_num);
}




/*
 * Return the minimum value of a device-allocated CudaGlobalPtr-array
 */

FLOAT thrustGetMinVal(CudaGlobalPtr<FLOAT> &diff2s)
{
	thrust::device_ptr<FLOAT> dp = thrust::device_pointer_cast(~diff2s);
	thrust::device_ptr<FLOAT> pos = thrust::min_element(dp, dp + diff2s.size);
	unsigned int pos_index = thrust::distance(dp, pos);
	FLOAT min_val;
	HANDLE_ERROR(cudaMemcpy(&min_val, &diff2s.d_ptr[pos_index], sizeof(FLOAT), cudaMemcpyDeviceToHost));
	return(min_val);
}


static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

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
inline
void mapWeights(CudaGlobalPtr<FLOAT> &mapped_weights, unsigned orientation_num, unsigned translation_num,
		HealpixSampling &sampling, long int ipart,
		std::vector< long unsigned > &iover_transes, std::vector< long unsigned > &ihiddens,
		std::vector< long unsigned > &iorientclasses, std::vector< long unsigned > &iover_rots,
		MultidimArray<FLOAT> &Mweight, unsigned long current_oversampling, unsigned long nr_trans)
{

	int nr_over_orient = sampling.oversamplingFactorOrientations(current_oversampling);
	int nr_over_trans =sampling.oversamplingFactorTranslations(current_oversampling);

	for (long unsigned i = 0; i < orientation_num; i++)
	{
		long unsigned iover_rot = iover_rots[i];
		for (long unsigned j = 0; j < translation_num; j++)
		{
			long unsigned iover_trans = iover_transes[j];
			long unsigned ihidden = iorientclasses[i] * nr_trans + ihiddens[j];
			long unsigned ihidden_over = ihidden * nr_over_orient * nr_over_trans + nr_over_trans * iover_rot + iover_trans;
			mapped_weights[(long unsigned) (i * translation_num + j)] =
					DIRECT_A2D_ELEM(Mweight, ipart, ihidden_over);
		}
	}
}

inline
long unsigned imageTranslation(
		CudaGlobalPtr<FLOAT> &Fimgs_real, CudaGlobalPtr<FLOAT> &Fimgs_imag,
		CudaGlobalPtr<FLOAT> &Fimgs_nomask_real, CudaGlobalPtr<FLOAT> &Fimgs_nomask_imag,
		long int itrans_min, long int itrans_max, int adaptive_oversampling , HealpixSampling &sampling,
		std::vector<double> &oversampled_translations_x, std::vector<double> &oversampled_translations_y, std::vector<double> &oversampled_translations_z,
		unsigned long nr_oversampled_trans, std::vector<MultidimArray<Complex> > &global_fftshifts_ab_current, std::vector<MultidimArray<Complex> > &global_fftshifts_ab2_current,
		MultidimArray<Complex > &local_Fimgs_shifted, MultidimArray<Complex > &local_Fimgs_shifted_nomask,
		std::vector< long unsigned > &iover_transes, std::vector< long unsigned > &itranses, std::vector< long unsigned > &ihiddens,
		unsigned image_size)
{

	long unsigned translation_num(0), ihidden(0);

	for (long int itrans = itrans_min, iitrans = 0; itrans <= itrans_max; itrans++, ihidden++)
	{
		sampling.getTranslations(itrans, adaptive_oversampling,
				oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);

		for (long int iover_trans = 0; iover_trans < nr_oversampled_trans; iover_trans++, iitrans++)
		{
			Complex* myAB;
			myAB = (adaptive_oversampling == 0 ) ? global_fftshifts_ab_current[iitrans].data : global_fftshifts_ab2_current[iitrans].data;


			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(local_Fimgs_shifted)
			{
				FLOAT a = (*(myAB + n)).real;
				FLOAT b = (*(myAB + n)).imag;

				// Fimg_shift
				FLOAT real = a * (DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted, n)).real
						- b *(DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted, n)).imag;
				FLOAT imag = a * (DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted, n)).imag
						+ b *(DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted, n)).real;
				Fimgs_real[translation_num * image_size + n] = real;
				Fimgs_imag[translation_num * image_size + n] = imag;

				// Fimg_shift_nomask
				real = a * (DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted_nomask, n)).real
						- b *(DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted_nomask, n)).imag;
				imag = a * (DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted_nomask, n)).imag
						+ b *(DIRECT_MULTIDIM_ELEM(local_Fimgs_shifted_nomask, n)).real;
				Fimgs_nomask_real[translation_num * image_size + n] = real;
				Fimgs_nomask_imag[translation_num * image_size + n] = imag;
			}

			translation_num ++;

			ihiddens.push_back(ihidden);
			itranses.push_back(itrans);
			iover_transes.push_back(iover_trans);
		}
	}

	Fimgs_real.size = translation_num * image_size;
	Fimgs_imag.size = translation_num * image_size;

	Fimgs_nomask_real.size = translation_num * image_size;
	Fimgs_nomask_imag.size = translation_num * image_size;

	return translation_num;
}


void generateEulerMatrices(
		FLOAT padding_factor,
		std::vector< double > &rots,
		std::vector< double > &tilts,
		std::vector< double > &psis,
		CudaGlobalPtr<FLOAT> &eulers,
		bool inverse)
{
	double alpha, beta, gamma;
    double ca, sa, cb, sb, cg, sg;
    double cc, cs, sc, ss;

	for (long int i = 0; i < rots.size(); i++)
	{
	    //TODO In a sense we're doing RAD2DEG just to do DEG2RAD here.
	    //The only place the degree value is actually used is in the metadata assignment.

	    alpha = DEG2RAD(rots[i]);
	    beta  = DEG2RAD(tilts[i]);
	    gamma = DEG2RAD(psis[i]);

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


long int generateProjectionSetup(
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		bool coarse,
		unsigned iclass,
		std::vector< double > &rots,
		std::vector< double > &tilts,
		std::vector< double > &psis,
		std::vector< long unsigned > &iorientclasses,
		std::vector< long unsigned > &iover_rots)
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
					iorientclasses.push_back(iorientclass);
					iover_rots.push_back(iover_rot);

					rots.push_back(oversampled_rot[iover_rot]);
					tilts.push_back(oversampled_tilt[iover_rot]);
					psis.push_back(oversampled_psi[iover_rot]);

					orientation_num ++;
				}
			}
		}
	}

	return orientation_num;
}


void runWavgKernel(
		Cuda3DProjectorKernel &projector,
		CudaGlobalPtr<FLOAT> &eulers,
		CudaGlobalPtr<FLOAT> &Fimgs_real,
	    CudaGlobalPtr<FLOAT> &Fimgs_imag,
	    CudaGlobalPtr<FLOAT> &Fimgs_nomask_real,
 	    CudaGlobalPtr<FLOAT> &Fimgs_nomask_imag,
 	    CudaGlobalPtr<FLOAT> &sorted_weights,
 	    CudaGlobalPtr<FLOAT> &ctfs,
 	    CudaGlobalPtr<FLOAT> &Minvsigma2s,
 	    CudaGlobalPtr<FLOAT> &wdiff2s_parts,
 	    CudaGlobalPtr<FLOAT> &wavgs_real,
	    CudaGlobalPtr<FLOAT> &wavgs_imag,
	    CudaGlobalPtr<FLOAT> &Fweights,
	    OptimisationParamters op,
	    MlOptimiser *baseMLO,
	    long unsigned orientation_num,
	    long unsigned translation_num,
	    unsigned image_size,
	    long int ipart,
	    int group_id,
	    int exp_iclass)
{
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	unsigned orientation_chunks = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);

	dim3 block_dim = splitCudaBlocks(orientation_chunks,false);

	CUDA_GPU_TIC("cuda_kernel_wavg");

	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);
	cuda_kernel_wavg<<<block_dim,BLOCK_SIZE>>>(
			~eulers,
			projector,
			image_size,
			orientation_num,
			~Fimgs_real, ~Fimgs_imag,
			~Fimgs_nomask_real, ~Fimgs_nomask_imag,
			~sorted_weights, ~ctfs, ~Minvsigma2s,
			~wdiff2s_parts,
			~wavgs_real,
			~wavgs_imag,
			~Fweights,
			translation_num,
			(FLOAT) op.sum_weight[ipart],
			(FLOAT) op.significant_weight[ipart],
			baseMLO->refs_are_ctf_corrected
			);

	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	float used = 100*((float)(total - avail)/(float)total);
	std::cerr << "Device memory used @ wavg: " << used << "%" << std::endl;
	CUDA_GPU_TAC("cuda_kernel_wavg");
}


void runDiff2KernelCoarse(
		Cuda3DProjectorKernel &projector,
		CudaGlobalPtr<FLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<FLOAT> &Fimgs_real,
		CudaGlobalPtr<FLOAT> &Fimgs_imag,
		CudaGlobalPtr<FLOAT> &eulers,
		CudaGlobalPtr<FLOAT> &diff2s,
		OptimisationParamters op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		unsigned image_size,
		int ipart,
		int group_id,
		int exp_iclass)
{

	CUDA_GPU_TIC("runProjAndDifferenceKernelCoarse");

	cuda_kernel_diff2_course<<<orientation_num,BLOCK_SIZE,translation_num*BLOCK_SIZE*sizeof(FLOAT)>>>(
			~eulers,
			~Fimgs_real,
			~Fimgs_imag,
			projector,
			~gpuMinvsigma2,
			~diff2s,
			translation_num,
			image_size,
			op.highres_Xi2_imgs[ipart] / 2.);


	CUDA_GPU_TAC("runProjAndDifferenceKernelCoarse");
}


void runDiff2KernelFine(
		Cuda3DProjectorKernel &projector,
		CudaGlobalPtr<FLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<FLOAT> &Fimgs_real,
		CudaGlobalPtr<FLOAT> &Fimgs_imag,
		CudaGlobalPtr<FLOAT> &eulers,
		CudaGlobalPtr<long unsigned> &rot_id,
		CudaGlobalPtr<long unsigned> &rot_idx,
		CudaGlobalPtr<long unsigned> &trans_idx,
		CudaGlobalPtr<long unsigned> &job_idx,
		CudaGlobalPtr<long unsigned> &job_num,
		CudaGlobalPtr<FLOAT> &diff2s,
		OptimisationParamters op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		long unsigned significant_num,
		unsigned image_size,
		int ipart,
		int group_id,
		int exp_iclass)
{
	CUDA_CPU_TIC("kernel_init_1");

	// Since we hijack Minvsigma to carry a bit more info into the GPU-kernel
	// we need to make a modified copy, since the global object shouldn't be
	// changed
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
	{
		gpuMinvsigma2[n] = *(op.local_Minvsigma2s[ipart].data + n );
	}

	if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
		{
			gpuMinvsigma2[n] *= (DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n));
		}
	}
	// TODO :    + Assure accuracy with the implemented GPU-based ctf-scaling
	//           + Make setting of myscale robust between here and above.
	//  (scale_correction turns off by default with only one group: ml_optimiser-line 1067,
	//   meaning small-scale test will probably not catch this malfunctioning when/if it breaks.)
	if (baseMLO->do_scale_correction)
	{
		FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
		{
			gpuMinvsigma2[n] *= (myscale*myscale);
		}
	}
    long int block_num = job_num.size;
    dim3 block_dim = splitCudaBlocks(block_num,false);

	CUDA_GPU_TIC("imagMemCp");
	gpuMinvsigma2.cp_to_device();
	Fimgs_real.put_on_device(translation_num * image_size);
	Fimgs_imag.put_on_device(translation_num * image_size);
	CUDA_GPU_TAC("imagMemCp");

	CUDA_GPU_TIC("pairListMemCp");
	rot_id.put_on_device(significant_num); //FIXME this is not used
	rot_idx.put_on_device(significant_num);
	trans_idx.put_on_device(significant_num);
	job_idx.put_on_device(block_num);
	job_num.put_on_device(block_num);
	CUDA_GPU_TAC("pairListMemCp");

	CUDA_CPU_TOC("kernel_init_1");

	CUDA_GPU_TIC("kernel_diff_proj");


// Could be used to automate __ldg() fallback runtime within cuda_kernel_diff2.
//				cudaDeviceProp dP;
//				cudaGetDeviceProperties(&dP, 0);
//				printf("-arch=sm_%d%d\n", dP.major, dP.minor);

	if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) // do cross-correlation instead of diff
	{
		// FIXME  make _CC
		printf("Cross correlation is not supported yet.");
		exit(0);
	}
	cuda_kernel_diff2_fine<<<block_dim,BLOCK_SIZE>>>(
			~eulers,
			~Fimgs_real,
			~Fimgs_imag,
			projector,
			~gpuMinvsigma2,
			~diff2s,
			image_size,
			op.highres_Xi2_imgs[ipart] / 2.,
			orientation_num,
			translation_num,
			block_num, //significant_num,
			~rot_idx,
			~trans_idx,
			~job_idx,
			~job_num);

	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	float used = 100*((float)(total - avail)/(float)total);
	std::cerr << "Device memory used @ diff2: " << used << "%" << std::endl;
	CUDA_GPU_TAC("kernel_diff_proj");
}
