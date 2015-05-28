#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_helper_kernels.cuh"
#include "src/gpu_utils/cuda_projection_kernels.cuh"
#include "src/gpu_utils/cuda_difference_kernels.cuh"
#include "src/gpu_utils/cuda_ProjDiff_kernels.cuh"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/complex.h"
#include <fstream>
#include <cuda.h>
#include "src/parallel.h"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>


static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

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

	for (long unsigned i = 0; i < orientation_num; i++)
	{
		long unsigned iover_rot = iover_rots[i];
		for (long unsigned j = 0; j < translation_num; j++)
		{
			long unsigned iover_trans = iover_transes[j];
			long unsigned ihidden = iorientclasses[i] * nr_trans + ihiddens[j];
			long unsigned ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden,
									  current_oversampling, iover_rot, iover_trans);
			mapped_weights[(long unsigned) i * translation_num + j] =
					DIRECT_A2D_ELEM(Mweight, ipart, ihidden_over);
			//Mweight[(i)*(v).xdim+(j)]
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
		    eulers[9 * i + 0] = ( cg * cc - sg * sa) * padding_factor; //00
		    eulers[9 * i + 1] = (-sg * cc - cg * sa) * padding_factor; //10
		    eulers[9 * i + 2] = ( sc )               * padding_factor; //20
		    eulers[9 * i + 3] = ( cg * cs + sg * ca) * padding_factor; //01
		    eulers[9 * i + 4] = (-sg * cs + cg * ca) * padding_factor; //11
		    eulers[9 * i + 5] = ( ss )               * padding_factor; //21
		    eulers[9 * i + 6] = (-cg * sb )          * padding_factor; //02
		    eulers[9 * i + 7] = ( sg * sb )          * padding_factor; //12
		    eulers[9 * i + 8] = ( cb )               * padding_factor; //22
		}
		else
		{
		    eulers[9 * i + 0] = ( cg * cc - sg * sa) * padding_factor; //00
		    eulers[9 * i + 1] = ( cg * cs + sg * ca) * padding_factor; //01
		    eulers[9 * i + 2] = (-cg * sb )          * padding_factor; //02
		    eulers[9 * i + 3] = (-sg * cc - cg * sa) * padding_factor; //10
		    eulers[9 * i + 4] = (-sg * cs + cg * ca) * padding_factor; //11
		    eulers[9 * i + 5] = ( sg * sb )          * padding_factor; //12
		    eulers[9 * i + 6] = ( sc )               * padding_factor; //20
		    eulers[9 * i + 7] = ( ss )               * padding_factor; //21
		    eulers[9 * i + 8] = ( cb )               * padding_factor; //22
		}
	}
}


void generateEulerMatrices(
		std::vector< double > &psis,
		CudaGlobalPtr<FLOAT> &eulers,
		bool inverse)
{
    double gamma, c, s;

	for (long int i = 0; i < psis.size(); i++)
	{
	    //TODO In a sense we're doing RAD2DEG just to do DEG2RAD here.
	    //The only place the degree value is actually used is in the metadata assignment.

	    gamma = DEG2RAD(psis[i]);
	    sincos(gamma, &s, &c);

		if(inverse) //Noticed here that inverse actually yields the opposite (Hmmm)
		{
		    eulers[4 * i + 0] =  c; //00
		    eulers[4 * i + 1] = -s; //10
		    eulers[4 * i + 3] =  s; //01
		    eulers[4 * i + 4] =  c; //11
		}
		else
		{
		    eulers[4 * i + 0] =  c; //00
		    eulers[4 * i + 1] =  s; //01
		    eulers[4 * i + 3] = -s; //10
		    eulers[4 * i + 4] =  c; //11
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

	unsigned parts_size(sp.nr_psi * sp.nr_oversampled_rot);
	std::vector< double > rots_parts(parts_size);
	std::vector< double > tilts_parts(parts_size);
	std::vector< double > psis_parts(parts_size);
	std::vector< long unsigned > iorientclasses_parts(parts_size);
	std::vector< long unsigned > iover_rots_parts(parts_size);

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
					iorientclasses_parts[ipart] = iorientclass;
					iover_rots_parts[ipart] = iover_rot;

					rots_parts[ipart] = oversampled_rot[iover_rot];
					tilts_parts[ipart] = oversampled_tilt[iover_rot];
					psis_parts[ipart] = oversampled_psi[iover_rot];

					orientation_num ++;
				}
			}
		}

		//TODO check that the following sort always works out

		if (sp.current_oversampling > 0)
		{
			int oversampling_per_psi = ROUND(std::pow(2., sp.current_oversampling));
			int oversampling_per_dir = ROUND(std::pow(4., sp.current_oversampling));

			//Sort the angles to have coalesced rot/tilt order
			for (unsigned i = 0; i < oversampling_per_dir; i++) //Loop over the perturbed dir pairs
			{
				for (unsigned j = 0; j < sp.nr_psi; j++)
				{
					for (unsigned k = 0; k < oversampling_per_psi; k++) //two psis per perturbed dir pair
					{
						unsigned ij = j*oversampling_per_psi*oversampling_per_dir + i*oversampling_per_psi + k;

						iorientclasses.push_back(iorientclasses_parts[ij]);
						iover_rots.push_back(iover_rots_parts[ij]);

						rots.push_back(rots_parts[ij]);
						tilts.push_back(tilts_parts[ij]);
						psis.push_back(psis_parts[ij]);
					}
				}
			}
		}
		else
		{
			for (unsigned i = 0; i < iorientclasses_parts.size(); i++)
			{
				iorientclasses.push_back(iorientclasses_parts[i]);
				iover_rots.push_back(iover_rots_parts[i]);

				rots.push_back(rots_parts[i]);
				tilts.push_back(tilts_parts[i]);
				psis.push_back(psis_parts[i]);
			}
		}
	}

	return orientation_num;
}

void generateModelProjections(
		CudaGlobalPtr<FLOAT > &model_real,
		CudaGlobalPtr<FLOAT > &model_imag,
		CudaGlobalPtr<FLOAT> &Frefs_real,
		CudaGlobalPtr<FLOAT> &Frefs_imag,
		CudaGlobalPtr<FLOAT> &eulers,
		long unsigned orientation_num,
		unsigned image_size,
		unsigned max_r,
		unsigned img_x,
		unsigned img_y,
		unsigned mdl_x,
		unsigned mdl_y,
		unsigned mdl_z,
		unsigned mdl_init_y,
		unsigned mdl_init_z)
{

	int max_r2 = max_r * max_r;
	int min_r2_nn = 0; // r_min_nn * r_min_nn;  //FIXME add nn-algorithm

	/*===========================
	 *      TEXTURE STUFF
	 * ==========================*/

	// create channel to describe data type (bits,bits,bits,bits,type)
	// TODO model should carry real & imag in separate channels of the same texture
	cudaChannelFormatDesc channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaArray*        modelArray_real;
	cudaArray* 		  modelArray_imag;
	cudaExtent        volumeSize = make_cudaExtent(mdl_x, mdl_y, mdl_z);

	//allocate device memory for cuda 3D array
	cudaMalloc3DArray(&modelArray_real, &channel, volumeSize);
	cudaMalloc3DArray(&modelArray_imag, &channel, volumeSize);

	//set cuda array copy parameters to be supplied to copy-command
	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent   = volumeSize;

	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_real;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_real.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);
	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_imag;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_imag.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);

	// Create texture object// Specify texture
    struct cudaResourceDesc resDesc_real,resDesc_imag;
    memset(&resDesc_real, 0, sizeof(resDesc_real));
    memset(&resDesc_imag, 0, sizeof(resDesc_imag));
    resDesc_real.resType = cudaResourceTypeArray;
    resDesc_imag.resType = cudaResourceTypeArray;
    resDesc_real.res.array.array = modelArray_real;
    resDesc_imag.res.array.array = modelArray_imag;

    struct cudaTextureDesc texDesc_real, texDesc_imag;
    memset(&texDesc_real, 0, sizeof(texDesc_real));
    memset(&texDesc_imag, 0, sizeof(texDesc_imag));
    for(int n=0; n<3; n++)
	{
    	texDesc_real.addressMode[n]=cudaAddressModeClamp;
    	texDesc_imag.addressMode[n]=cudaAddressModeClamp;
	}
    texDesc_real.filterMode       = cudaFilterModeLinear;
    texDesc_real.readMode         = cudaReadModeElementType;
    texDesc_real.normalizedCoords = false;
    texDesc_imag.filterMode       = cudaFilterModeLinear;
    texDesc_imag.readMode         = cudaReadModeElementType;
    texDesc_real.normalizedCoords = false;

	cudaTextureObject_t texModel_real = 0;
	cudaCreateTextureObject(&texModel_real, &resDesc_real, &texDesc_real, NULL);
	cudaTextureObject_t texModel_imag = 0;
	cudaCreateTextureObject(&texModel_imag, &resDesc_imag, &texDesc_imag, NULL);

	Frefs_real.size = orientation_num * image_size;
	Frefs_real.device_alloc();
	Frefs_imag.size = orientation_num * image_size;
	Frefs_imag.device_alloc();

	unsigned int orient1, orient2;
	if(orientation_num>65535)
	{
		orient1 = ceil(sqrt(orientation_num));
		orient2 = orient1;
	}
	else
	{
		orient1 = orientation_num;
		orient2 = 1;
	}

	dim3 block_dim(orient1,orient2);
	std::cerr << "using block dimensions " << orient1 << "," << orient2 <<  std::endl;

#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)
	// we CAN use read-associated interpolation (fast, inaccurate)...
	cuda_kernel_PAV_TTI<<<block_dim,BLOCK_SIZE>>>(
													~eulers,
													~Frefs_real,
													~Frefs_imag,
													texModel_real,
													texModel_imag,
													max_r,
													max_r2,
													min_r2_nn,
													image_size,
													orientation_num,
													img_x,
													img_y,
													mdl_init_y,
													mdl_init_z);
	cudaDestroyTextureObject(texModel_real);
	cudaDestroyTextureObject(texModel_imag);
	cudaFreeArray(modelArray_real);
	cudaFreeArray(modelArray_imag);
#elif !defined(CUDA_DOUBLE_PRECISION)	// ...or explicit interpolation (slow, accurate)
	cuda_kernel_PAV_TTE<<<block_dim,BLOCK_SIZE>>>(
													~eulers,
													~Frefs_real,
													~Frefs_imag,
													texModel_real,
													texModel_imag,
													max_r,
													max_r2,
													min_r2_nn,
													image_size,
													orientation_num,
													img_x,
													img_y,
													mdl_init_y,
													mdl_init_z);

	cudaDestroyTextureObject(texModel_real);
	cudaDestroyTextureObject(texModel_imag);
	cudaFreeArray(modelArray_real);
	cudaFreeArray(modelArray_imag);
#else // under double precision, texture won't work.
	model_real.device_alloc();
    model_real.cp_to_device();
 	model_imag.device_alloc();
    model_imag.cp_to_device();
	cuda_kernel_PAV_TGE<<<block_dim,BLOCK_SIZE>>>(
													~model_real,
													~model_imag,
													~eulers,
													~Frefs_real,
													~Frefs_imag,
													max_r,
													max_r2,
													min_r2_nn,
													image_size,
													orientation_num,
													img_x,
													img_y,
													mdl_x,
													mdl_y,
													mdl_init_y,
													mdl_init_z);
	model_real.free_device();
	model_imag.free_device();
#endif

	//unbind texture reference to free resource



}


void runDifferenceKernel(CudaGlobalPtr<FLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<FLOAT > &Fimgs_real,
		CudaGlobalPtr<FLOAT > &Fimgs_imag,
		CudaGlobalPtr<FLOAT > &Frefs_real,
		CudaGlobalPtr<FLOAT > &Frefs_imag,
		CudaGlobalPtr<long unsigned > &rotidx,
		CudaGlobalPtr<long unsigned > &transidx,
		CudaGlobalPtr<long unsigned > &ihidden_overs,
		OptimisationParamters &op,
		MlOptimiser *baseMLO,
		long unsigned translation_num,
		long unsigned orientation_num,
		long unsigned significant_num,
		unsigned image_size,
		int ipart,
		int group_id,
		CudaGlobalPtr<FLOAT > &diff2s
		)
{
	/*====================================
	   Initiate Particle Related On GPU
	======================================*/
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

	gpuMinvsigma2.cp_to_device();

	Fimgs_real.size = translation_num * image_size;
	Fimgs_real.device_alloc();
	Fimgs_real.cp_to_device();
	Fimgs_imag.size = translation_num * image_size;
	Fimgs_imag.device_alloc();
	Fimgs_imag.cp_to_device();
	rotidx.size = significant_num;
	rotidx.device_alloc();
	rotidx.cp_to_device();
	transidx.size = significant_num;
	transidx.device_alloc();
	transidx.cp_to_device();
	ihidden_overs.size = significant_num;
	ihidden_overs.device_alloc();
	ihidden_overs.cp_to_device();

	/*====================================
				Kernel Calls
	======================================*/
	unsigned orient1, orient2;

	if(significant_num>65535)
	{
		orient1 = ceil(sqrt(significant_num));
		orient2 = orient1;
	}
	else
	{
		orient1 = significant_num;
		orient2 = 1;
	}
	dim3 block_dim(orient1,orient2);

	CUDA_CPU_TOC("kernel_init_1");
	CUDA_GPU_TIC("kernel_diff_noproj");
	// Could be used to automate __ldg() fallback runtime within cuda_kernel_diff2.
//				cudaDeviceProp dP;
//				cudaGetDeviceProperties(&dP, 0);
//				printf("-arch=sm_%d%d\n", dP.major, dP.minor);

	if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) // do cross-correlation instead of diff
	{
		cuda_kernel_D2_CC<<<block_dim,BLOCK_SIZE>>>(~Frefs_real, ~Frefs_imag, ~Fimgs_real, ~Fimgs_imag, ~gpuMinvsigma2,  ~diff2s,
														image_size, op.highres_Xi2_imgs[ipart],
														significant_num,
														translation_num,
														~rotidx,
														~transidx);
	}
	else
	{
		cuda_kernel_D2<<<block_dim,BLOCK_SIZE>>>(~Frefs_real, ~Frefs_imag, ~Fimgs_real, ~Fimgs_imag, ~gpuMinvsigma2, ~diff2s,
													image_size, op.highres_Xi2_imgs[ipart] / 2.,
													significant_num,
													translation_num,
													~rotidx,
													~transidx,
													~ihidden_overs);
	}
	CUDA_GPU_TAC("kernel_diff_noproj");
	HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here
	CUDA_GPU_TOC("kernel_diff_noproj");
	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	float used = 100*((float)(total - avail)/(float)total);
	std::cerr << "Device memory used @ diff2: " << used << "%" << std::endl;
}

#if !defined(CUDA_DOUBLE_PRECISION)
void runProjAndDifferenceKernel(
		CudaGlobalPtr<FLOAT > &model_real,
		CudaGlobalPtr<FLOAT > &model_imag,
		CudaGlobalPtr<FLOAT > &gpuMinvsigma2,
		CudaGlobalPtr<FLOAT> &Fimgs_real,
		CudaGlobalPtr<FLOAT> &Fimgs_imag,
		CudaGlobalPtr<FLOAT> &eulers,
		CudaGlobalPtr<long unsigned> &rotidx,
		CudaGlobalPtr<long unsigned> &transidx,
		CudaGlobalPtr<long unsigned> &trans_num,
		CudaGlobalPtr<long unsigned> &ihidden_overs,
		CudaGlobalPtr<FLOAT> &diff2s,
		OptimisationParamters op,
		MlOptimiser *baseMLO,
		long unsigned orientation_num,
		long unsigned translation_num,
		long unsigned significant_num,
		unsigned image_size,
		unsigned max_r,
		int ipart,
		int group_id,
		int exp_iclass)
{

	CUDA_CPU_TIC("kernel_init_1");
	int max_r2 = max_r * max_r;
	int min_r2_nn = 0; // r_min_nn * r_min_nn;  //FIXME add nn-algorithm

	/*===========================
	 *      TEXTURE STUFF
	 * ==========================*/

	// create channel to describe data type (bits,bits,bits,bits,type)
	// TODO model should carry real & imag in separate channels of the same texture
	cudaChannelFormatDesc channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaArray*        modelArray_real;
	cudaArray* 		  modelArray_imag;
	cudaExtent        volumeSize = make_cudaExtent(baseMLO->mymodel.PPref[exp_iclass].data.xdim,
			   	   	   	   	   	   	   	   	       baseMLO->mymodel.PPref[exp_iclass].data.ydim,
			   	   	   	   	   	   	   	   	       baseMLO->mymodel.PPref[exp_iclass].data.zdim);

	//allocate device memory for cuda 3D array
	cudaMalloc3DArray(&modelArray_real, &channel, volumeSize);
	cudaMalloc3DArray(&modelArray_imag, &channel, volumeSize);

	//set cuda array copy parameters to be supplied to copy-command
	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent   = volumeSize;

	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_real;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_real.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);
	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_imag;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_imag.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);

	// Create texture object// Specify texture
	struct cudaResourceDesc resDesc_real,resDesc_imag;
	memset(&resDesc_real, 0, sizeof(resDesc_real));
	memset(&resDesc_imag, 0, sizeof(resDesc_imag));
	resDesc_real.resType = cudaResourceTypeArray;
	resDesc_imag.resType = cudaResourceTypeArray;
	resDesc_real.res.array.array = modelArray_real;
	resDesc_imag.res.array.array = modelArray_imag;

	struct cudaTextureDesc texDesc_real, texDesc_imag;
	memset(&texDesc_real, 0, sizeof(texDesc_real));
	memset(&texDesc_imag, 0, sizeof(texDesc_imag));
	for(int n=0; n<3; n++)
	{
		texDesc_real.addressMode[n]=cudaAddressModeClamp;
		texDesc_imag.addressMode[n]=cudaAddressModeClamp;
	}
	texDesc_real.filterMode       = cudaFilterModeLinear;
	texDesc_real.readMode         = cudaReadModeElementType;
	texDesc_real.normalizedCoords = false;
	texDesc_imag.filterMode       = cudaFilterModeLinear;
	texDesc_imag.readMode         = cudaReadModeElementType;
	texDesc_real.normalizedCoords = false;

	cudaTextureObject_t texModel_real = 0;
	cudaCreateTextureObject(&texModel_real, &resDesc_real, &texDesc_real, NULL);
	cudaTextureObject_t texModel_imag = 0;
	cudaCreateTextureObject(&texModel_imag, &resDesc_imag, &texDesc_imag, NULL);


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

	gpuMinvsigma2.cp_to_device();

	unsigned orient1, orient2;
    unsigned long block_num = trans_num.size;//significant_num;

	if(block_num>65535)
	{
		orient1 = ceil(sqrt(block_num));
		orient2 = orient1;
	}
	else
	{
		orient1 = block_num;
		orient2 = 1;
	}
	dim3 block_dim(orient1,orient2);

	Fimgs_real.size = translation_num * image_size;
	Fimgs_real.device_alloc();
	Fimgs_real.cp_to_device();
	Fimgs_imag.size = translation_num * image_size;
	Fimgs_imag.device_alloc();
	Fimgs_imag.cp_to_device();
	rotidx.size = block_num;
	rotidx.device_alloc();
	rotidx.cp_to_device();
	transidx.size = block_num;
	transidx.device_alloc();
	transidx.cp_to_device();
//  trans_num.size = block_num; // already set
    trans_num.device_alloc();
    trans_num.cp_to_device();
	ihidden_overs.size = block_num;
	ihidden_overs.device_alloc();
	ihidden_overs.cp_to_device();


	CUDA_CPU_TOC("kernel_init_1");
	CUDA_GPU_TIC("kernel_diff_proj");


// Could be used to automate __ldg() fallback runtime within cuda_kernel_diff2.
//				cudaDeviceProp dP;
//				cudaGetDeviceProperties(&dP, 0);
//				printf("-arch=sm_%d%d\n", dP.major, dP.minor);

	if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) // do cross-correlation instead of diff
	{
		// FIXME  make _CC
		exit(0);
//		cuda_kernel_PAV_TTI_D2_CC<<<block_dim,BLOCK_SIZE>>>(~eulers,
//														 ~Fimgs_real,
//														 ~Fimgs_imag,
//														 texModel_real,
//													 	 texModel_imag,
//														 ~gpuMinvsigma2,
//														 ~diff2s,
//														 image_size,
//														 op.highres_Xi2_imgs[ipart] / 2.,
//														 orientation_num,
//														 translation_num,
//														 significant_num,
//														 ~rotidx,
//														 ~transidx,
//		 	 	 	 	 	 	 	 	 	 	 	 	 ~trans_num,
//														 ~ihidden_overs,
//														 max_r,
//													     max_r2,
//													     min_r2_nn,
//														 op.local_Minvsigma2s[0].xdim,
//														 op.local_Minvsigma2s[0].ydim,
//														 baseMLO->mymodel.PPref[exp_iclass].data.yinit,
//														 baseMLO->mymodel.PPref[exp_iclass].data.zinit);
//		cudaDestroyTextureObject(texModel_real);
//		cudaDestroyTextureObject(texModel_imag);
//		cudaFreeArray(modelArray_real);
//		cudaFreeArray(modelArray_imag);
	}
	else
	{
		cuda_kernel_PAV_TTI_D2<<<block_dim,BLOCK_SIZE>>>(~eulers,
														 ~Fimgs_real,
														 ~Fimgs_imag,
														 texModel_real,
														 texModel_imag,
														 ~gpuMinvsigma2,
														 ~diff2s,
														 image_size,
														 op.highres_Xi2_imgs[ipart] / 2.,
														 orientation_num,
														 translation_num,
														 block_num, //significant_num,
														 ~rotidx,
														 ~transidx,
														 ~trans_num,
														 ~ihidden_overs,
														 max_r,
														 max_r2,
														 min_r2_nn,
														 op.local_Minvsigma2s[0].xdim,
														 op.local_Minvsigma2s[0].ydim,
														 baseMLO->mymodel.PPref[exp_iclass].data.yinit,
														 baseMLO->mymodel.PPref[exp_iclass].data.zinit);
		size_t avail;
		size_t total;
		cudaMemGetInfo( &avail, &total );
		float used = 100*((float)(total - avail)/(float)total);
		std::cerr << "Device memory used @ diff2: " << used << "%" << std::endl;
		cudaDestroyTextureObject(texModel_real);
		cudaDestroyTextureObject(texModel_imag);
		cudaFreeArray(modelArray_real);
		cudaFreeArray(modelArray_imag);
	}
	CUDA_GPU_TAC("kernel_diff_proj");
	HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here

	CUDA_GPU_TOC("kernel_diff_proj");

}
#endif


#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_FETCH_COUNT 4

__global__ void cuda_kernel_backproject(
		int *g_xs,
		int *g_ys,
		int *g_zs,
		FLOAT *g_model_real,
		FLOAT *g_model_imag,
		FLOAT *g_weight,
		FLOAT *g_eulers,
		FLOAT *g_wavgs_real,
		FLOAT *g_wavgs_imag,
		FLOAT *g_Fweights,
		int max_r2, FLOAT scale2,
		unsigned img_xy, unsigned long img_count, unsigned img_x, unsigned img_y,
		unsigned mdl_x, unsigned mdl_y, int mdl_inity, int mdl_initz,
		int N)
{
	unsigned gid = threadIdx.x / 4;
	unsigned mid = threadIdx.x % 4;
	unsigned gm = gid * 4 + mid;
	unsigned pit = (gid * 4 + mid)*BACKPROJECTION4_FETCH_COUNT;
	unsigned global_idx = blockIdx.x * BACKPROJECTION4_GROUP_SIZE + gid;

	int X(0),Y(0),Z(0);

	if (global_idx < N)
	{
		X = g_xs[global_idx];
		Y = g_ys[global_idx];
		Z = g_zs[global_idx];
	}
	else
		X = mdl_x * 10; // Padding coordinate, place outside images

	int ax(0), ay(0);

	if (mid == 1)
		ax = 1;
	else if (mid == 2)
		ay = 1;
	else if (mid == 3)
	{
		ax = 1;
		ay = 1;
	}

	bool  is_neg_x;
	FLOAT d, w;
	FLOAT xp,yp,zp;
	int x,y,idx;

	__shared__ FLOAT s_e[BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT];

	__shared__ FLOAT s_weight[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ FLOAT s_value_real[BACKPROJECTION4_GROUP_SIZE*4];
	__shared__ FLOAT s_value_imag[BACKPROJECTION4_GROUP_SIZE*4];

	s_weight[gm] = 0.0f;
	s_value_real[gm] = 0.0f;
	s_value_imag[gm] = 0.0f;

	for (int img = 0, b = BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT; img < img_count; img ++, b += 9)
	{
		if (b+9 > BACKPROJECTION4_BLOCK_SIZE*BACKPROJECTION4_FETCH_COUNT)
		{
			__syncthreads();

			int img_9 = img*9+pit;
			if (img_9 < img_count*9)
			{
				s_e[pit+0] = g_eulers[img_9+0];
				s_e[pit+1] = g_eulers[img_9+1];
				s_e[pit+2] = g_eulers[img_9+2];
				s_e[pit+3] = g_eulers[img_9+3];
			}

			__syncthreads();
			b = 0;
		}

		zp = (s_e[b+6] * X + s_e[b+7] * Y + s_e[b+8] * Z);

		if (fabsf(zp) > 0.87f) continue; //Within the unit cube, sqrt(3)/2=0.866

		yp = (s_e[b+3] * X + s_e[b+4] * Y + s_e[b+5] * Z);
		xp = (s_e[b+0] * X + s_e[b+1] * Y + s_e[b+2] * Z);

		if (xp < 0.0f)
		{
			yp = -yp;
			xp = -xp;
			is_neg_x = true;
		}
		else
			is_neg_x = false;

		x = (int) floorf(xp) + ax;
		y = (int) floorf(yp) + ay;

		if (x * x + y * y > max_r2) continue;

		if (y < 0 && x == 0)
		{
			is_neg_x = !is_neg_x;
			y = -y;
		}

		xp = (s_e[b+0] * x + s_e[b+3] * y) * scale2;
		yp = (s_e[b+1] * x + s_e[b+4] * y) * scale2;
		zp = (s_e[b+2] * x + s_e[b+5] * y) * scale2;

		if (xp < 0.0f) //Flip sign
		{
			xp = fabsf(X+xp);
			yp = fabsf(Y+yp);
			zp = fabsf(Z+zp);
		}
		else
		{
			xp = fabsf(X-xp);
			yp = fabsf(Y-yp);
			zp = fabsf(Z-zp);
		}

		if (xp < 1.0f && yp < 1.0f && zp < 1.0f)
		{
			if (y < 0) y += img_y;
			idx = img*img_xy + y * img_x + x;
			w = g_Fweights[idx];

			if (w > 0.0f)
			{
				d = (1.0f - xp) * (1.0f - yp) * (1.0f - zp);

				s_weight[gm] += w * d;
				s_value_real[gm] += g_wavgs_real[idx] * d;
				if (is_neg_x) s_value_imag[gm] -= g_wavgs_imag[idx] * d;
				else          s_value_imag[gm] += g_wavgs_imag[idx] * d;
			}
		}
	}

	__syncthreads();

	if (mid == 0)
	{
		FLOAT sum = s_weight[gid*4 + 0] + s_weight[gid*4 + 1] + s_weight[gid*4 + 2] + s_weight[gid*4 + 3];
		if (sum != 0.0f)
			g_weight[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] = sum;
	}
	else if (mid == 1)
	{
		FLOAT sum = s_value_real[gid*4 + 0] + s_value_real[gid*4 + 1] + s_value_real[gid*4 + 2] + s_value_real[gid*4 + 3];
		if (sum != 0.0f)
			g_model_real[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] = sum;
	}
	else if (mid == 2)
	{
		FLOAT sum = s_value_imag[gid*4 + 0] + s_value_imag[gid*4 + 1] + s_value_imag[gid*4 + 2] + s_value_imag[gid*4 + 3];
		if (sum != 0.0f)
			g_model_imag[(Z-mdl_initz)*mdl_x*mdl_y + (Y-mdl_inity)*mdl_x + X] = sum;
	}
}

static void runBackprojectKernel(
		CudaGlobalPtr<FLOAT> &wavgs_real,
		CudaGlobalPtr<FLOAT> &wavgs_imag,
		CudaGlobalPtr<FLOAT> &Fweights,
		CudaGlobalPtr<FLOAT> &eulers,
		CudaGlobalPtr<FLOAT> &model_real,
		CudaGlobalPtr<FLOAT> &model_imag,
		CudaGlobalPtr<FLOAT> &weight,
		int max_r, FLOAT scale2, //grid scale 2D -> 3D squared
		int img_xy, long img_count, int img_x, int img_y,
		int mdl_x, int mdl_y, int mdl_z, int mdl_inity, int mdl_initz)
{
	int max_r2 = max_r * max_r;

	CudaGlobalPtr<int> xs(mdl_x*mdl_y*mdl_z); // >52% will actually be used, allocate some padding
	CudaGlobalPtr<int> ys(xs.size);
	CudaGlobalPtr<int> zs(xs.size);
	unsigned N(0);

	for (int x = 0; x < mdl_x; x ++)
	{
		for (int y = mdl_inity; y < mdl_y; y++)
		{
			for (int z = mdl_initz; z < mdl_z; z++)
			{
				if (x*x + y*y + z*z <= max_r2 * scale2 * 1.2f)
				{
					xs[N] = x;
					ys[N] = y;
					zs[N] = z;
					N ++;
				}
			}
		}
	}
	xs.size = N + N%BACKPROJECTION4_GROUP_SIZE;
	ys.size = xs.size;
	zs.size = xs.size;

	xs.device_alloc();
	ys.device_alloc();
	zs.device_alloc();

	xs.cp_to_device();
	ys.cp_to_device();
	zs.cp_to_device();

	int grid_dim = ceil((float)N / BACKPROJECTION4_GROUP_SIZE);
	dim3 block_dim( BACKPROJECTION4_GROUP_SIZE *4 );

	cuda_kernel_backproject<<<grid_dim,block_dim>>>(
			~xs,~ys,~zs,
			~model_real,
			~model_imag,
			~weight,
			~eulers,
			~wavgs_real,
			~wavgs_imag,
			~Fweights,
			max_r2,
			scale2,
			img_xy,
			img_count,
			img_x,
			img_y,
			mdl_x,
			mdl_y,
			mdl_inity,
			mdl_initz,
			N);
}









void MlOptimiserCuda::doThreadExpectationSomeParticles(unsigned thread_id)
{
	size_t first_ipart = 0, last_ipart = 0;
	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			unsigned my_ori_particle = baseMLO->exp_my_first_ori_particle + ipart;
			SamplingParameters sp;
			sp.nr_particles = baseMLO->mydata.ori_particles[my_ori_particle].particles_id.size();

			OptimisationParamters op(sp.nr_particles, my_ori_particle);

			// In the first iteration, multiple seeds will be generated
			// A single random class is selected for each pool of images, and one does not marginalise over the orientations
			// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
			// If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

			// Decide which classes to integrate over (for random class assignment in 1st iteration)
			sp.iclass_max = baseMLO->mymodel.nr_classes - 1;
			// low-pass filter again and generate the seeds
			if (baseMLO->do_generate_seeds)
			{
				if (baseMLO->do_firstiter_cc && baseMLO->iter == 1)
				{
					// In first (CC) iter, use a single reference (and CC)
					sp.iclass_min = sp.iclass_max = 0;
				}
				else if ( (baseMLO->do_firstiter_cc && baseMLO->iter == 2) ||
						(!baseMLO->do_firstiter_cc && baseMLO->iter == 1))
				{
					// In second CC iter, or first iter without CC: generate the seeds
					// Now select a single random class
					// exp_part_id is already in randomized order (controlled by -seed)
					// WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
					sp.iclass_min = sp.iclass_max = divide_equally_which_group(baseMLO->mydata.numberOfOriginalParticles(), baseMLO->mymodel.nr_classes, op.my_ori_particle);
				}
			}
			// Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
			for (long int iori = baseMLO->exp_my_first_ori_particle; iori <= baseMLO->exp_my_last_ori_particle; iori++)
			{
				if (iori == my_ori_particle) break;
				op.metadata_offset += baseMLO->mydata.ori_particles[iori].particles_id.size();
			}

			baseMLO->getFourierTransformsAndCtfs(my_ori_particle, op.metadata_offset, op.Fimgs, op.Fimgs_nomask, op.Fctfs,
					op.old_offset, op.prior, op.power_imgs, op.highres_Xi2_imgs,
					op.pointer_dir_nonzeroprior, op.pointer_psi_nonzeroprior, op.directions_prior, op.psi_prior);

			if (baseMLO->do_realign_movies && baseMLO->movie_frame_running_avg_side > 0)
			{
				baseMLO->calculateRunningAveragesOfMovieFrames(my_ori_particle, op.Fimgs, op.power_imgs, op.highres_Xi2_imgs);
			}

			// To deal with skipped alignments/rotations
			if (baseMLO->do_skip_align)
			{
				sp.itrans_min = sp.itrans_max = sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
						my_ori_particle - baseMLO->exp_my_first_ori_particle;
			}
			else
			{
				sp.itrans_min = 0;
				sp.itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;

				if (baseMLO->do_skip_rotate)
				{
					sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
							my_ori_particle - baseMLO->exp_my_first_ori_particle;
				}
				else
				{
					sp.idir_min = sp.ipsi_min = 0;
					sp.idir_max = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) - 1;
					sp.ipsi_max = baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior ) - 1;
				}
			}

			// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
			op.significant_weight.resize(sp.nr_particles, -1.);

			// Only perform a second pass when using adaptive oversampling
			int nr_sampling_passes = (baseMLO->adaptive_oversampling > 0) ? 2 : 1;

			for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
			{
				if (baseMLO->strict_highres_exp > 0.)
					// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
					sp.current_image_size = baseMLO->coarse_size;
				else if (baseMLO->adaptive_oversampling > 0)
					// Use smaller images in the first pass, larger ones in the second pass
					sp.current_image_size = (ipass == 0) ? baseMLO->coarse_size : baseMLO->mymodel.current_size;
				else
					sp.current_image_size = baseMLO->mymodel.current_size;

				// Use coarse sampling in the first pass, oversampled one the second pass
				sp.current_oversampling = (ipass == 0) ? 0 : baseMLO->adaptive_oversampling;

				sp.nr_dir = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior);
				sp.nr_psi = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
				sp.nr_trans = (baseMLO->do_skip_align) ? 1 : baseMLO->sampling.NrTranslationalSamplings();
				sp.nr_oversampled_rot = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				sp.nr_oversampled_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);

				CUDA_CPU_TIC("getAllSquaredDifferences");
				getAllSquaredDifferences(ipass, op, sp);
				CUDA_CPU_TOC("getAllSquaredDifferences");
				CUDA_CPU_TIC("convertAllSquaredDifferencesToWeights");
				convertAllSquaredDifferencesToWeights(ipass, op, sp);
				CUDA_CPU_TOC("convertAllSquaredDifferencesToWeights");
			}

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;

			CUDA_CPU_TIC("storeWeightedSums");
			storeWeightedSums(op, sp);
			CUDA_CPU_TOC("storeWeightedSums");
		}
	}
}





void MlOptimiserCuda::getAllSquaredDifferences(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp)
{

	CUDA_CPU_TIC("diff_pre_gpu");

	//for scale_correction
	int group_id;

	//printf("sp.nr_oversampled_rot=%d\n", (unsigned)sp.nr_oversampled_rot);

	op.Mweight.resize(sp.nr_particles, baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);
	op.Mweight.initConstant(-999.);
	if (exp_ipass==0)
	{
		op.Mcoarse_significant.clear();
	}

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	MultidimArray<Complex > Fref;
	Fref.resize(op.local_Minvsigma2s[0]);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{
		if (baseMLO->mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			CudaGlobalPtr<FLOAT> gpuMinvsigma2(image_size);
			gpuMinvsigma2.device_alloc();

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			std::vector< double > rots, tilts, psis;

			CUDA_CPU_TIC("projection_1");
			CUDA_CPU_TIC("generateProjectionSetup");
			long unsigned orientation_num = generateProjectionSetup(
					op,
					sp,
					baseMLO,
					exp_ipass == 0, //coarse
					exp_iclass,
					rots, tilts, psis,
					iorientclasses,
					iover_rots);

			CUDA_CPU_TOC("generateProjectionSetup");
			CUDA_CPU_TIC("generateEulerMatrices");
			CudaGlobalPtr<FLOAT> eulers(9 * orientation_num);

			generateEulerMatrices(
					baseMLO->mymodel.PPref[exp_iclass].padding_factor,
					rots,
					tilts,
					psis,
					eulers,
					!IS_NOT_INV);

		    eulers.device_alloc();
			eulers.cp_to_device();
			CUDA_CPU_TOC("generateEulerMatrices");
			CUDA_CPU_TIC("modelAssignment");

//			CudaGlobalPtr<FLOAT > model_real;
//			float* holder = new float[(baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim];
//			model_real.h_ptr = holder;
//			model_real.size=(baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim;
//			model_real.h_do_free = true;

			CudaGlobalPtr<FLOAT > model_real((baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim);
			CudaGlobalPtr<FLOAT > model_imag((baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim);

			for(unsigned i = 0; i < model_real.size; i++)
			{
				model_real[i] = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].real;
				model_imag[i] = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].imag;
			}

			CudaGlobalPtr<FLOAT> Frefs_real;
		    CudaGlobalPtr<FLOAT> Frefs_imag;

			CUDA_CPU_TOC("modelAssignment");
			bool do_combineProjAndDiff = true; //TODO add control flag
			if(!do_combineProjAndDiff)
			{
				CUDA_CPU_TIC("generateModelProjections_diff");
				generateModelProjections(
						model_real,
						model_imag,
						Frefs_real,
						Frefs_imag,
						eulers,
						orientation_num,
						image_size,
						XMIPP_MIN(baseMLO->mymodel.PPref[exp_iclass].r_max, op.local_Minvsigma2s[0].xdim - 1),
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						baseMLO->mymodel.PPref[exp_iclass].data.xdim,
						baseMLO->mymodel.PPref[exp_iclass].data.ydim,
						baseMLO->mymodel.PPref[exp_iclass].data.zdim,
						baseMLO->mymodel.PPref[exp_iclass].data.yinit,
						baseMLO->mymodel.PPref[exp_iclass].data.zinit);
				model_real.free_device();
				model_imag.free_device();
				eulers.free_device();
				CUDA_CPU_TOC("generateModelProjections_diff");
			}
			CUDA_CPU_TOC("projection_1");

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CUDA_CPU_TIC("translation_1");

				CudaGlobalPtr<FLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans);
				CudaGlobalPtr<FLOAT> Fimgs_imag(image_size * sp.nr_trans * sp.nr_oversampled_trans);

				long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					baseMLO->sampling.getTranslations(itrans, sp.current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (sp.current_oversampling == 0)
						{
							myAB = (Fref.ydim == baseMLO->coarse_size) ? baseMLO->global_fftshifts_ab_coarse[itrans].data
									: baseMLO->global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * sp.nr_oversampled_trans +  iover_trans;
							myAB = (baseMLO->strict_highres_exp > 0.) ? baseMLO->global_fftshifts_ab2_coarse[iitrans].data
									: baseMLO->global_fftshifts_ab2_current[iitrans].data;
						}


						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
						{
							FLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
							FLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (baseMLO->do_scale_correction)
							{
								//group_id = mydata.getGroupId(part_id);
								FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
							}
							Fimgs_real[translation_num * image_size + n] = real;
							Fimgs_imag[translation_num * image_size + n] = imag;
						}
						translation_num ++;

						ihiddens.push_back(ihidden);
						itranses.push_back(itrans);
						iover_transes.push_back(iover_trans);
					}
				}

				CUDA_CPU_TOC("translation_1");

				/*===========================================
				   Determine significant comparison indices
				=============================================*/
				//      This section is annoying to test because
				//		it can't complete on first pass, since
				//		the significance has never been set


				CUDA_CPU_TIC("pair_list_1");

				CudaGlobalPtr<long unsigned> transidx(orientation_num*translation_num), rotidx(orientation_num*translation_num);
				CudaGlobalPtr<long unsigned> ihidden_overs(orientation_num*translation_num);
				CudaGlobalPtr<long unsigned> trans_num(orientation_num*translation_num);
				long unsigned coarse_num = sp.nr_dir*sp.nr_psi*sp.nr_trans;
				long unsigned significant_num(0);
//				long int check_num=0;
				long unsigned k=0;
				if (exp_ipass == 0)
				{
					op.Mcoarse_significant.resize(coarse_num, 1);
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						trans_num[k]=0;
						transidx[k]=translation_num+1;//set higher than max(j) so that XMIPP_MIN() sets
						int tk=0;
						for (long unsigned j = 0; j < translation_num; j++)
						{
							ihidden_overs[significant_num] = i * sp.nr_trans + j;
							if(tk>=PROJDIFF_CHUNK_SIZE)
							{
								tk=0;             // reset counter
//								check_num+=trans_num[k];
								k++;              // use new element
								trans_num[k]=0;   // prepare next element for ++ incrementing
								transidx[k]=translation_num+1; //set higher than max(j) so that XMIPP_MIN() sets
							}
							tk++;                 // increment limit
							trans_num[k]++;       // increment number of transes this ProjDiff-block
							rotidx[k] = i;
							transidx[k] = XMIPP_MIN(j,transidx[k]);
							significant_num++;
						}
//						check_num+=trans_num[k];
						k++;   // use new element
					}
					trans_num.size=k;
				}
				else
				{
					for (long unsigned i = 0; i < orientation_num; i++)
					{
						trans_num[k]=0;
						transidx[k]=translation_num+1;//set higher than max(j) so that XMIPP_MIN() sets
						int tk=0;

						long int iover_rot = iover_rots[i];
						long int coarse_rot = floor(i/sp.nr_oversampled_rot);
						for (long unsigned j = 0; j < translation_num; j++)
						{
							long int iover_trans = iover_transes[j];
							long int coarse_trans = floor(j/sp.nr_oversampled_trans);
							long int ihidden = iorientclasses[i] * sp.nr_trans + ihiddens[j];

							if(DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden)==1)
							{
								ihidden_overs[significant_num] = baseMLO->sampling.getPositionOversampledSamplingPoint(ihidden,
										                  sp.current_oversampling, iover_rot, iover_trans);
								if(tk>=PROJDIFF_CHUNK_SIZE)
								{
									tk=0;             // reset counter
//									check_num+=trans_num[k];
									k++;              // use new element
									trans_num[k]=0;   // prepare next element for ++ incrementing
									transidx[k]=translation_num+1; //set higher than max(j) so that XMIPP_MIN() sets
								}
								tk++;                 // increment limit
								trans_num[k]++;       // increment number of transes this ProjDiff-block
								rotidx[k] = i;
								transidx[k] = XMIPP_MIN(j,transidx[k]);
								significant_num++;
							}
							else if(tk!=0) // start a new one - we expect transes to be sequential.
							{
								tk=0;             // reset counter
//								check_num+=trans_num[k];
								k++;              // use new element
								trans_num[k]=0;   // prepare next element for ++ incrementing
								transidx[k]=translation_num+1; //set higher than max(j) so that XMIPP_MIN() sets
							}

						}
//						check_num+=trans_num[k];
						k++;   // use new element
					}
					trans_num.size=k;
				}
				//  check_num should equal significant_num here, and be less or equal to  PROJDIFF_CHUNK_SIZE*trans_num.size

				CUDA_CPU_TOC("pair_list_1");
//				std::cerr << "orientation_num "<< orientation_num << std::endl;
//				std::cerr << "translation_num "<< translation_num << std::endl;
//				std::cerr << "my_nr_significant_coarse_samples "<< DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NR_SIGN) << std::endl;
//				std::cerr << "significant_num "<< significant_num << std::endl;

				CudaGlobalPtr<FLOAT> diff2s(orientation_num*translation_num);
				diff2s.device_alloc();

#if !defined(CUDA_DOUBLE_PRECISION)
				if(do_combineProjAndDiff)
				{
					runProjAndDifferenceKernel(model_real,
											   model_imag,
											   gpuMinvsigma2,
										       Fimgs_real,
										       Fimgs_imag,
										       eulers,
										       rotidx,
										       transidx,
										       trans_num,
										       ihidden_overs,
										       diff2s,
										       op,
										       baseMLO,
										       orientation_num,
										       translation_num,
										       significant_num,
										       image_size,
											    XMIPP_MIN(baseMLO->mymodel.PPref[exp_iclass].r_max, op.local_Minvsigma2s[0].xdim - 1),
										       ipart,
										       group_id,
										       exp_iclass
											 );
					eulers.free_device();
				}
				else
#endif
				{
					runDifferenceKernel(gpuMinvsigma2,
										Fimgs_real,
										Fimgs_imag,
										Frefs_real,
										Frefs_imag,
										rotidx,
										transidx,
										ihidden_overs,
										op,
										baseMLO,
										translation_num,
										orientation_num,
										significant_num,
										image_size,
										ipart,
										group_id,
										diff2s
										);
				}
				/*====================================
				    	   Retrieve Results
				======================================*/

				diff2s.cp_to_host(); // FIXME may not be needed since we copy it back in ConvetToWeights()
//				for (long unsigned k = 0; k < 100; k++)
//				{
//					std::cerr << diff2s[k] << std::endl;
//				}
				if (exp_ipass == 0)
				{
					op.Mcoarse_significant.clear();
				}

				/*====================================
				    	Write To Destination
				======================================*/


				CUDA_CPU_TIC("collect_data_1");

//				freopen(text,"w",stdout);
				long unsigned m=0;
				for (long unsigned k = 0; k < trans_num.size; k++)
				{
					for (int itrans=0;  itrans < trans_num[k]; itrans++, m++)
					{
						long unsigned i = rotidx[k];
						long unsigned j = transidx[k]+itrans;
						double diff2 = diff2s[i * translation_num + j];
//						printf("%4.8f \n",DIRECT_A2D_ELEM(op.Mweight, ipart, n));
//						printf("%4.8f, %i, %i \n",diff2,i,j);// << std::endl;

						DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_overs[m]) = diff2; // TODO if we can write diff2 to the correct pos in the kernel we can just memcpy to a pointer and use thrust to find min
						// Keep track of minimum of all diff2, only for the last image in this series
						if (diff2 < op.min_diff2[ipart])
							op.min_diff2[ipart] = diff2;
					}
				}
//				fclose(stdout);

				CUDA_CPU_TOC("collect_data_1");
			} // end loop ipart
		} // end if class significant
	} // end loop iclass
}

void MlOptimiserCuda::convertAllSquaredDifferencesToWeights(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp)
{
	op.sum_weight.clear();
	op.sum_weight.resize(sp.nr_particles, 0.);

	// loop over all particles inside this ori_particle
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		double exp_thisparticle_sumweight = 0.;

		double old_offset_z;
		double old_offset_x = XX(op.old_offset[ipart]);
		double old_offset_y = YY(op.old_offset[ipart]);
		if (baseMLO->mymodel.data_dim == 3)
			old_offset_z = ZZ(op.old_offset[ipart]);

		if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		{
			// Binarize the squared differences array to skip marginalisation
			double mymindiff2 = 99.e10;
			long int myminidx = -1;
			// Find the smallest element in this row of op.Mweight
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
			{

				double cc = DIRECT_A2D_ELEM(op.Mweight, ipart, i);
				// ignore non-determined cc
				if (cc == -999.)
					continue;

				// just search for the maximum
				if (cc < mymindiff2)
				{
					mymindiff2 = cc;
					myminidx = i;
				}
			}
			// Set all except for the best hidden variable to zero and the smallest element to 1
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
				DIRECT_A2D_ELEM(op.Mweight, ipart, i)= 0.;

			DIRECT_A2D_ELEM(op.Mweight, ipart, myminidx)= 1.;
			exp_thisparticle_sumweight += 1.;

		}
		else
		{
			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{

				// Make PdfOffset calculation much faster...
				double myprior_x, myprior_y, myprior_z;
				if (baseMLO->mymodel.ref_dim == 2)
				{
					myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(op.prior[ipart]);
					myprior_y = YY(op.prior[ipart]);
					if (baseMLO->mymodel.data_dim == 3)
						myprior_z = ZZ(op.prior[ipart]);
				}

				/*=========================================
						Fetch+generate Orientation data
				===========================================*/
				CudaGlobalPtr<FLOAT >  pdf_orientation(sp.nr_dir * sp.nr_psi);
				pdf_orientation.size = sp.nr_dir * sp.nr_psi;
				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						//std::cerr << "orient "  << idir << "," << iorient <<  std::endl;
						// Get prior for this direction
						if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
						{
							pdf_orientation[iorient] = baseMLO->mymodel.pdf_class[exp_iclass];
						}
						else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						{
							pdf_orientation[iorient] = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							// P(orientation) = P(idir|dir_prior) * P(ipsi|psi_prior)
							// This is the probability of the orientation, given the gathered
							// statistics of all assigned orientations of the dataset, since we
							// are assigning a gaussian prior to all parameters.
							pdf_orientation[iorient] = op.directions_prior[idir] * op.psi_prior[ipsi];
						}
					}
				}
//				long int ihidden = iorientclass * sp.nr_trans;

				/*=========================================
						Fetch+generate Translation data
				===========================================*/
				CudaGlobalPtr<FLOAT >  pdf_offset(sp.nr_trans);

				int jtrans=0;
				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++,jtrans++)
				{
					//std::cerr << "trans " << itrans << "," << jtrans <<  std::endl;
			        // To speed things up, only calculate pdf_offset at the coarse sampling.
					// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
					double offset_x = old_offset_x + baseMLO->sampling.translations_x[itrans];
					double offset_y = old_offset_y + baseMLO->sampling.translations_y[itrans];
					double tdiff2 = (offset_x - myprior_x) * (offset_x - myprior_x) + (offset_y - myprior_y) * (offset_y - myprior_y);
					if (baseMLO->mymodel.data_dim == 3)
					{
						double offset_z = old_offset_z + baseMLO->sampling.translations_z[itrans];
						tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
					}
					// P(offset|sigma2_offset)
					// This is the probability of the offset, given the model offset and variance.
					if (baseMLO->mymodel.sigma2_offset < 0.0001)
						pdf_offset[jtrans] = ( tdiff2 > 0.) ? 0. : 1.;
					else
						pdf_offset[jtrans] = exp ( tdiff2 / (-2. * baseMLO->mymodel.sigma2_offset) ) / ( 2. * PI * baseMLO->mymodel.sigma2_offset );
				}

// TODO : Put back when  convertAllSquaredDifferencesToWeights is GPU-parallel.
//							// TMP DEBUGGING
//							if (baseMLO->mymodel.orientational_prior_mode != NOPRIOR && (pdf_offset==0. || pdf_orientation==0.))
//							{
//								pthread_mutex_lock(&global_mutex);
//								std::cerr << " pdf_offset= " << pdf_offset << " pdf_orientation= " << pdf_orientation << std::endl;
//								std::cerr << " ipart= " << ipart << " part_id= " << part_id << std::endl;
//								std::cerr << " iorient= " << iorient << " idir= " << idir << " ipsi= " << ipsi << std::endl;
//								//std::cerr << " sp.nr_psi= " << sp.nr_psi << " exp_nr_dir= " << exp_nr_dir << " sp.nr_trans= " << sp.nr_trans << std::endl;
//								for (long int i = 0; i < op.directions_prior.size(); i++)
//									std::cerr << " op.directions_prior["<<i<<"]= " << op.directions_prior[i] << std::endl;
//								for (long int i = 0; i < op.psi_prior.size(); i++)
//									std::cerr << " op.psi_prior["<<i<<"]= " << op.psi_prior[i] << std::endl;
//								REPORT_ERROR("ERROR! pdf_offset==0.|| pdf_orientation==0.");
//								//pthread_mutex_unlock(&global_mutex);
//							}
//							if (sp.nr_oversampled_rot == 0)
//								REPORT_ERROR("sp.nr_oversampled_rot == 0");
//							if (sp.nr_oversampled_trans == 0)
//								REPORT_ERROR("sp.nr_oversampled_trans == 0");

				// Now first loop over iover_rot, because that is the order in op.Mweight as well
//				long int ihidden_over = ihidden * sp.nr_oversampled_rot * sp.nr_oversampled_trans;

				/*=========================================
					  Kernel call over all combinations
				===========================================*/

				// One block will be started for each (coarse) orientation, and will process all (coarse) transes,
				// and since oversmapling is by factors of 2 on 5 dofs, we get 2^5=32 fine comparisons per coarse.
				// In case of higher oversampling this is simply factors of 32, making it warp-perfect. Having 21
				// coarse transes allows a block to finish in 21, 100% utilized, warp passes.

				int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;

				bool do_gpu_sumweight = true;  //TODO add control flag
				if(oversamples>=SUM_BLOCK_SIZE && do_gpu_sumweight) // Send task to GPU where warps can access automatically coalesced oversamples
				{
					//std::cerr << "summing weights on GPU... baseMLO->mymodel.pdf_class[exp_iclass] = " << baseMLO->mymodel.pdf_class[sp.iclass_min] <<  std::endl;
					pdf_orientation.device_alloc();
					pdf_orientation.cp_to_device();
					pdf_offset.device_alloc();
					pdf_offset.cp_to_device();

					CudaGlobalPtr<FLOAT >  thisparticle_sumweight(sp.nr_dir * sp.nr_psi);  // This will be reduced in a second step.
					thisparticle_sumweight.device_alloc();

					CudaGlobalPtr<FLOAT >  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim]),
													sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples);
					Mweight.device_alloc();
					Mweight.cp_to_device();

					dim3 block_dim(sp.nr_dir,sp.nr_psi);
					//std::cerr << "using block dimensions " << sp.nr_dir << "," << sp.nr_psi <<  std::endl;
					cuda_kernel_sumweight_oversampling<<<block_dim,SUM_BLOCK_SIZE>>>(	~pdf_orientation,
																						~pdf_offset,
																						~Mweight,
																						~thisparticle_sumweight,
																						op.min_diff2[ipart],
																						sp.nr_trans,
																						oversamples
																					 );

					Mweight.cp_to_host(); //FIXME make wider in scope; pass to storeWsums() to be used in collect-step. Needs som coordination with else() below.
					Mweight.free_device();  //FIXME see line above
					thisparticle_sumweight.cp_to_host();
					thisparticle_sumweight.free_device();

					// The reduced entity *MUST* be double to avoid loss of information// TODO better reduction
					for (long int n = 0; n < sp.nr_dir * sp.nr_psi; n++)
					{
						exp_thisparticle_sumweight += (double)thisparticle_sumweight[n];
					}
				}
				else // Not enough oversamples to utilize GPU resources effciently with current CUDA-kernel.
				{
					//std::cerr << "summing weights on CPU... " <<  std::endl;
					for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
					{
						for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
						{
							long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;
							long int ihidden = iorientclass * sp.nr_trans;
							for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
							{
								long int ihidden_over = ihidden * sp.nr_oversampled_rot * sp.nr_oversampled_trans;
								for (long int iover_rot = 0; iover_rot < sp.nr_oversampled_rot; iover_rot++)
								{
									// Then loop over iover_trans
									for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, ihidden_over++)
									{
										// Only exponentiate for determined values of op.Mweight
										// (this is always true in the first pass, but not so in the second pass)
										// Only deal with this sampling point if its weight was significant
										if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) < 0.)
										{
											DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) = 0.;
										}
										else
										{
											// Set the weight base to the probability of the parameters given the prior
											double weight = pdf_orientation[iorient] * pdf_offset[itrans];
											double diff2 = DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) - op.min_diff2[ipart];
											// next line because of numerical precision of exp-function
											if (diff2 > 700.) weight = 0.;
											// TODO: use tabulated exp function?
											else weight *= exp(-diff2);
											// Store the weight
											DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_over) = weight;

											// Keep track of sum and maximum of all weights for this particle
											// Later add all to exp_thisparticle_sumweight, but inside this loop sum to local thisthread_sumweight first
											exp_thisparticle_sumweight += weight;
										} // end if/else op.Mweight < 0.
									} // end loop iover_trans
								}// end loop iover_rot
							} // end loop itrans
						} // end loop ipsi
					} // end loop idir
				}                            //endif do_gpu_sumweight
			} // end loop exp_iclass
		} // end if iter==1

		//Store parameters for this particle
		op.sum_weight[ipart] = exp_thisparticle_sumweight;

#if defined(DEBUG_CUDA) && defined(__linux__)
		if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
		{
			printf("DEBUG_ERROR: zero sum of weights.\n");
			exit( EXIT_FAILURE );
		}
#endif

	} // end loop ipart

	if (exp_ipass==0)
	{
		op.Mcoarse_significant.resize(sp.nr_particles, XSIZE(op.Mweight));
	}

	CUDA_CPU_TIC("convert_post_kernel");
	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of op.sum_weight
	op.significant_weight.clear();
	op.significant_weight.resize(sp.nr_particles, 0.);
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		MultidimArray<FLOAT> sorted_weight;
		// Get the relevant row for this particle
		op.Mweight.getRow(ipart, sorted_weight);

		// Only select non-zero probabilities to speed up sorting
		long int np = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
		{
			if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
			{
				DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
				np++;
			}
		}
		sorted_weight.resize(np);

		// Sort from low to high values
		CUDA_CPU_TIC("sort");
#if defined(USE_THRUST) // Thrust seems incredibly slow in debug build this is clearly a FIXME
		thrust::sort(sorted_weight.data, sorted_weight.data + np);
#else
		sorted_weight.sort();
#endif
		CUDA_CPU_TOC("sort");

		double frac_weight = 0.;
		double my_significant_weight;
		long int my_nr_significant_coarse_samples = 0;
		for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
		{
			if (exp_ipass==0) my_nr_significant_coarse_samples++;
			my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
			frac_weight += my_significant_weight;
			if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
				break;
		}

		if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
		{
			std::cerr << " ipart= " << ipart << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
			std::cerr << " frac-weight= " << frac_weight << std::endl;
			std::cerr << " op.sum_weight[ipart]= " << op.sum_weight[ipart] << std::endl;
			Image<FLOAT> It;
			std::cerr << " XSIZE(op.Mweight)= " << XSIZE(op.Mweight) << std::endl;
			It()=op.Mweight;
			It() *= 10000;
			It.write("Mweight2.spi");
			std::cerr << "written Mweight2.spi" << std::endl;
			std::cerr << " np= " << np << std::endl;
			It()=sorted_weight;
			It() *= 10000;
			std::cerr << " XSIZE(sorted_weight)= " << XSIZE(sorted_weight) << std::endl;
			if (XSIZE(sorted_weight) > 0)
			{
				It.write("sorted_weight.spi");
				std::cerr << "written sorted_weight.spi" << std::endl;
			}
			REPORT_ERROR("my_nr_significant_coarse_samples == 0");
		}

		if (exp_ipass==0)
		{
			// Store nr_significant_coarse_samples for this particle
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NR_SIGN) = (double)my_nr_significant_coarse_samples;

			// Keep track of which coarse samplings were significant were significant for this particle
			for (int ihidden = 0; ihidden < XSIZE(op.Mcoarse_significant); ihidden++)
			{
				if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden) >= my_significant_weight)
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = true;
				else
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = false;
			}

		}
		op.significant_weight[ipart] = my_significant_weight;
	} // end loop ipart
	CUDA_CPU_TOC("convert_post_kernel");

}

// __global__ void cuda_kernel_wavg_fast   // REMOVED in commit
#if !defined(CUDA_DOUBLE_PRECISION)
__global__ void cuda_kernel_ProjAndWavg(
		FLOAT *g_eulers,
		cudaTextureObject_t texModel_real,
		cudaTextureObject_t texModel_imag,
		unsigned my_r_max,
		int max_r2,
		int min_r2_nn,
		unsigned image_size,
		unsigned long orientation_num,
	 	long int XSIZE_img,
	 	long int YSIZE_img,
	 	long int STARTINGY_mdl,
	 	long int STARTINGZ_mdl,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		FLOAT *g_imgs_nomask_real,
		FLOAT *g_imgs_nomask_imag,
		FLOAT* g_weights,
		FLOAT* g_ctfs,
		FLOAT* g_Minvsigma2s,
		FLOAT *g_wdiff2s_parts,
		FLOAT *g_wavgs_real,
		FLOAT *g_wavgs_imag,
		FLOAT* g_Fweights,
		unsigned long translation_num,
		FLOAT weight_norm,
		FLOAT significant_weight,
		bool refs_are_ctf_corrected)
{
	FLOAT xp, yp, zp;
	long int r2;
	bool is_neg_x;
	FLOAT ref_real, ref_imag;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
//	if( bid < orientation_num )
//	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  )),pixel;
		FLOAT Fweight;
		__shared__ FLOAT s_wavgs_real[BLOCK_SIZE];
		__shared__ FLOAT s_wavgs_imag[BLOCK_SIZE];
		__shared__ FLOAT s_wdiff2s_parts[BLOCK_SIZE];
		__shared__ FLOAT s_Minvsigma2s[BLOCK_SIZE];
		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			s_wavgs_real[tid]  = 0.0f;
			s_wavgs_imag[tid]  = 0.0f;
			s_wdiff2s_parts[tid] = 0.0f;
			Fweight = 0;

			pixel = pass * BLOCK_SIZE + tid;
			s_Minvsigma2s[tid]=g_Minvsigma2s[pixel];

			if(pixel<image_size)
			{
				unsigned long ref_pixel = bid * image_size + pixel;
				// Now istead of loading pre-calculated ref, we project it out from the texture-model
				//----------------------------------------------------------------------------------- =>
				int x = pixel % XSIZE_img;
				int y = (int)floorf( (float)pixel / (float)XSIZE_img);

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= YSIZE_img - my_r_max)
						y = y - YSIZE_img ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					xp = __ldg(&g_eulers[bid*9])   * x + __ldg(&g_eulers[bid*9+1]) * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = __ldg(&g_eulers[bid*9+3]) * x + __ldg(&g_eulers[bid*9+4]) * y;  // compared to CPU-based projection. This
					zp = __ldg(&g_eulers[bid*9+6]) * x + __ldg(&g_eulers[bid*9+7]) * y;  // propagates to dx00, dx10, and so on.
					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}
					yp -= STARTINGY_mdl;
					zp -= STARTINGZ_mdl;

					ref_real=tex3D<FLOAT>(texModel_real,xp+0.5f,yp+0.5f,zp+0.5f);
					ref_imag=tex3D<FLOAT>(texModel_imag,xp+0.5f,yp+0.5f,zp+0.5f);


					if (is_neg_x)
					{
						ref_imag = -ref_imag;
					}

				}
				else
				{
					ref_real=0.0f;
					ref_imag=0.0f;
				}
				//-----------------------------------------------------------------------------------  <=
				if (refs_are_ctf_corrected) //FIXME Create two kernels for the different cases
				{
					ref_real *= __ldg(&g_ctfs[pixel]);
					ref_imag *= __ldg(&g_ctfs[pixel]);
				}

				for (unsigned long itrans = 0; itrans < translation_num; itrans++)
				{
					FLOAT weight = __ldg(&g_weights[bid * translation_num + itrans]);

					if (weight >= significant_weight)
					{
						weight /= weight_norm;

						unsigned long img_pixel_idx = itrans * image_size + pixel;

						FLOAT diff_real = ref_real - g_imgs_real[img_pixel_idx];    // TODO  Put in texture (in such a way that fetching of next image might hit in cache)
						FLOAT diff_imag = ref_imag - g_imgs_imag[img_pixel_idx];

						s_wdiff2s_parts[tid] += weight * (diff_real*diff_real + diff_imag*diff_imag);

						FLOAT weightxinvsigma2 = weight * __ldg(&g_ctfs[pixel]) * s_Minvsigma2s[tid];

						s_wavgs_real[tid] += g_imgs_nomask_real[img_pixel_idx] * weightxinvsigma2;    // TODO  Put in texture (in such a way that fetching of next image might hit in cache)
						s_wavgs_imag[tid] += g_imgs_nomask_imag[img_pixel_idx] * weightxinvsigma2;

						Fweight += weightxinvsigma2 * __ldg(&g_ctfs[pixel]);
					}
				}
				g_wavgs_real[ref_pixel] += s_wavgs_real[tid];
				g_wavgs_imag[ref_pixel] += s_wavgs_imag[tid];
				g_wdiff2s_parts[ref_pixel] = s_wdiff2s_parts[tid]; //TODO this could be further reduced in here
				g_Fweights[ref_pixel] += Fweight; //TODO should be buffered into shared
			}
		}
//	}
}
#endif

void runWavgKernel(CudaGlobalPtr<FLOAT> &Frefs_real,
				   CudaGlobalPtr<FLOAT> &Frefs_imag,
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

	unsigned orient1, orient2;
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	unsigned orientation_chunks = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);
	if(orientation_chunks>65535)
	{
		orient1 = ceil(sqrt(orientation_chunks));
		orient2 = orient1;
	}
	else
	{
		orient1 = orientation_chunks;
		orient2 = 1;
	}
	dim3 block_dim(orient1,orient2);

	CUDA_GPU_TIC("cuda_kernel_wavg");

	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);
	cuda_kernel_wavg<<<block_dim,BLOCK_SIZE>>>(
										~Frefs_real, ~Frefs_imag, ~Fimgs_real, ~Fimgs_imag,
										~Fimgs_nomask_real, ~Fimgs_nomask_imag,
										~sorted_weights, ~ctfs, ~Minvsigma2s,
										~wdiff2s_parts,
										~wavgs_real,
										~wavgs_imag,
										~Fweights,
										orientation_num,
										translation_num,
										(FLOAT) op.sum_weight[ipart],
										(FLOAT) op.significant_weight[ipart],
										image_size,
										baseMLO->refs_are_ctf_corrected
										);
	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	float used = 100*((float)(total - avail)/(float)total);
	std::cerr << "Device memory used @ wavg: " << used << "%" << std::endl;
	CUDA_GPU_TAC("cuda_kernel_wavg");

	HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here

	CUDA_GPU_TOC("cuda_kernel_wavg");

//	Fimgs_real.free_device();
//	Fimgs_imag.free_device();
	Fimgs_nomask_real.free_device();
	Fimgs_nomask_imag.free_device();

	sorted_weights.free_device();
	ctfs.free_device();
	Minvsigma2s.free_device();
}


#if !defined(CUDA_DOUBLE_PRECISION)
void runProjAndWavgKernel(
		CudaGlobalPtr<FLOAT> &model_real,
		CudaGlobalPtr<FLOAT> &model_imag,
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
	    unsigned max_r,
	    long int ipart,
	    int group_id,
	    int exp_iclass
	   )
{
	int max_r2 = max_r * max_r;
	int min_r2_nn = 0; // r_min_nn * r_min_nn;  //FIXME add nn-algorithm

	/*===========================
	 *      TEXTURE STUFF
	 * ==========================*/

	// create channel to describe data type (bits,bits,bits,bits,type)
	// TODO model should carry real & imag in separate channels of the same texture
	cudaChannelFormatDesc channel = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
	cudaArray*        modelArray_real;
	cudaArray* 		  modelArray_imag;
	cudaExtent        volumeSize = make_cudaExtent(baseMLO->mymodel.PPref[exp_iclass].data.xdim,
			   	   	   	   	   	   	   	   	       baseMLO->mymodel.PPref[exp_iclass].data.ydim,
			   	   	   	   	   	   	   	   	       baseMLO->mymodel.PPref[exp_iclass].data.zdim);

	//allocate device memory for cuda 3D array
	cudaMalloc3DArray(&modelArray_real, &channel, volumeSize);
	cudaMalloc3DArray(&modelArray_imag, &channel, volumeSize);

	//set cuda array copy parameters to be supplied to copy-command
	cudaMemcpy3DParms copyParams = {0};
	copyParams.extent   = volumeSize;

	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_real;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_real.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);
	copyParams.kind     = cudaMemcpyHostToDevice;
	copyParams.dstArray = modelArray_imag;
	copyParams.srcPtr   = make_cudaPitchedPtr(model_imag.h_ptr,volumeSize.width*sizeof(FLOAT), volumeSize.height, volumeSize.depth);
	cudaMemcpy3D(&copyParams);

	// Create texture object// Specify texture
    struct cudaResourceDesc resDesc_real,resDesc_imag;
    memset(&resDesc_real, 0, sizeof(resDesc_real));
    memset(&resDesc_imag, 0, sizeof(resDesc_imag));
    resDesc_real.resType = cudaResourceTypeArray;
    resDesc_imag.resType = cudaResourceTypeArray;
    resDesc_real.res.array.array = modelArray_real;
    resDesc_imag.res.array.array = modelArray_imag;

    struct cudaTextureDesc texDesc_real, texDesc_imag;
    memset(&texDesc_real, 0, sizeof(texDesc_real));
    memset(&texDesc_imag, 0, sizeof(texDesc_imag));
    for(int n=0; n<3; n++)
	{
    	texDesc_real.addressMode[n]=cudaAddressModeClamp;
    	texDesc_imag.addressMode[n]=cudaAddressModeClamp;
	}
    texDesc_real.filterMode       = cudaFilterModeLinear;
    texDesc_real.readMode         = cudaReadModeElementType;
    texDesc_real.normalizedCoords = false;
    texDesc_imag.filterMode       = cudaFilterModeLinear;
    texDesc_imag.readMode         = cudaReadModeElementType;
    texDesc_real.normalizedCoords = false;

	cudaTextureObject_t texModel_real = 0;
	cudaCreateTextureObject(&texModel_real, &resDesc_real, &texDesc_real, NULL);
	cudaTextureObject_t texModel_imag = 0;
	cudaCreateTextureObject(&texModel_imag, &resDesc_imag, &texDesc_imag, NULL);
	unsigned orient1, orient2;
	//We only want as many blocks as there are chunks of orientations to be treated
	//within the same block (this is done to reduce memory loads in the kernel).
	unsigned orientation_chunks = orientation_num;//ceil((float)orientation_num/(float)REF_GROUP_SIZE);
	if(orientation_chunks>65535)
	{
		orient1 = ceil(sqrt(orientation_chunks));
		orient2 = orient1;
	}
	else
	{
		orient1 = orientation_chunks;
		orient2 = 1;
	}
	dim3 block_dim(orient1,orient2);

	CUDA_GPU_TIC("cuda_kernel_wavg");

	//cudaFuncSetCacheConfig(cuda_kernel_wavg_fast, cudaFuncCachePreferShared);
	cuda_kernel_ProjAndWavg<<<block_dim,BLOCK_SIZE>>>(~eulers,
													  texModel_real,
													  texModel_imag,
													  max_r,
													  max_r2,
													  min_r2_nn,
													  image_size,
													  orientation_num,
													  op.local_Minvsigma2s[0].xdim,
													  op.local_Minvsigma2s[0].ydim,
													  baseMLO->mymodel.PPref[exp_iclass].data.yinit,
													  baseMLO->mymodel.PPref[exp_iclass].data.zinit,
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
	cudaDestroyTextureObject(texModel_real);
	cudaDestroyTextureObject(texModel_imag);
	cudaFreeArray(modelArray_real);
	cudaFreeArray(modelArray_imag);

	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	float used = 100*((float)(total - avail)/(float)total);
	std::cerr << "Device memory used @ wavg: " << used << "%" << std::endl;
	CUDA_GPU_TAC("cuda_kernel_wavg");

	HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Apparently this is not required here

	CUDA_GPU_TOC("cuda_kernel_wavg");

//	Fimgs_real.free_device();
//	Fimgs_imag.free_device();
	Fimgs_nomask_real.free_device();
	Fimgs_nomask_imag.free_device();

	sorted_weights.free_device();
	ctfs.free_device();
	Minvsigma2s.free_device();
}
#endif

void MlOptimiserCuda::storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp)
{
	CUDA_CPU_TIC("store_pre_gpu");

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(true, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, op.Fimgs_nomask, op.Fctfs, op.local_Fimgs_shifted, op.local_Fimgs_shifted_nomask,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the op.local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);
		DIRECT_MULTIDIM_ELEM(op.local_Minvsigma2s[ipart], 0) = 1. / (baseMLO->sigma2_fudge * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], 0));
	}

	// Initialise the maximum of all weights to a negative value
	op.max_weight.clear();
	op.max_weight.resize(sp.nr_particles, -1.);

	// For norm_correction and scale_correction of all particles of this ori_particle
	std::vector<double> exp_wsum_norm_correction;
	std::vector<MultidimArray<double> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<MultidimArray<double> > thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(sp.nr_particles, 0.);

	// For scale_correction
	if (baseMLO->do_scale_correction)
	{
		MultidimArray<double> aux;
		aux.initZeros(baseMLO->mymodel.ori_size/2 + 1);
		exp_wsum_scale_correction_XA.resize(sp.nr_particles, aux);
		exp_wsum_scale_correction_AA.resize(sp.nr_particles, aux);
		thr_wsum_signal_product_spectra.resize(baseMLO->mymodel.nr_groups, aux);
		thr_wsum_reference_power_spectra.resize(baseMLO->mymodel.nr_groups, aux);
	}


	std::vector<double> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	Matrix2D<double> A;
	MultidimArray<FLOAT > Fimg_real, Fimg_imag;
	MultidimArray<Complex > Fimg, Fimg_otfshift_nomask;  //TODO remove, currently needed for Fourier stuff, which is based on the complex class
	MultidimArray<double> Fweight, Minvsigma2, Mctf;
	bool have_warned_small_scale = false;

	Fimg_real.resize(op.Fimgs[0]);
	Fimg_imag.resize(op.Fimgs[0]);
	Fimg.resize(op.Fimgs[0]);
	Fweight.resize(op.Fimgs[0]);

	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(op.Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(op.Fimgs[0]);
	Minvsigma2.initConstant(1.);

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_pdf_direction;
	std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	double thr_wsum_sigma2_offset;
	MultidimArray<double> thr_metadata, zeroArray;
	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(baseMLO->mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(baseMLO->mymodel.nr_groups, zeroArray);
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(baseMLO->sampling.NrDirections());
	thr_wsum_pdf_direction.resize(baseMLO->mymodel.nr_classes, zeroArray);
	// sumw_group is a double for each group
	thr_sumw_group.resize(baseMLO->mymodel.nr_groups, 0.);
	// wsum_pdf_class is a double for each class
	thr_wsum_pdf_class.resize(baseMLO->mymodel.nr_classes, 0.);
	if (baseMLO->mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(baseMLO->mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(baseMLO->mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a double
	thr_wsum_sigma2_offset = 0.;

	unsigned image_size = op.Fimgs[0].nzyxdim;

	CUDA_CPU_TOC("store_pre_gpu");

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{

		/*=======================================================================================
		                            REFERENCE PROJECTION GENERATION
		=======================================================================================*/

		// Since we will need the euler matrices for all projections in the data_collect stage,
		// we might as well make it wider in scope and retain it on the GPU until then. When we
		// switch from pair to bool, there won't be any need to remake it every class, but for
		// now we create only those matrices corresponding to significant orientations, which IS  * class-specific *

		std::vector< long unsigned > iorientclasses, iover_rots;
		std::vector< double > rots, tilts, psis;

		CUDA_CPU_TIC("projection_2");

		long unsigned orientation_num = generateProjectionSetup(
					op,
					sp,
					baseMLO,
					false,  //coarse
					exp_iclass,
					rots, tilts, psis,
					iorientclasses,
					iover_rots);


		CudaGlobalPtr<FLOAT> eulers(9 * orientation_num);

		generateEulerMatrices(
				baseMLO->mymodel.PPref[exp_iclass].padding_factor,
				rots,
				tilts,
				psis,
				eulers,
				!IS_NOT_INV);

	    eulers.device_alloc();
		eulers.cp_to_device();

		CudaGlobalPtr<FLOAT > model_real((baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim);
		CudaGlobalPtr<FLOAT > model_imag((baseMLO->mymodel.PPref[exp_iclass]).data.nzyxdim);

		for(unsigned i = 0; i < model_real.size; i++)
		{
			model_real[i] = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].real;
			model_imag[i] = (FLOAT) baseMLO->mymodel.PPref[exp_iclass].data.data[i].imag;
		}

		CudaGlobalPtr<FLOAT> Frefs_real;
		CudaGlobalPtr<FLOAT> Frefs_imag;

		bool do_combineProjAndWavg = true; //TODO add control flag
#if !defined(CUDA_DOUBLE_PRECISION)
		if(!do_combineProjAndWavg)
#endif
		{
			CUDA_CPU_TIC("generateModelProjections_wavg");
			generateModelProjections(
					model_real,
					model_imag,
					Frefs_real,
					Frefs_imag,
					eulers,
					orientation_num,
					image_size,
					XMIPP_MIN(baseMLO->mymodel.PPref[exp_iclass].r_max, op.local_Minvsigma2s[0].xdim - 1),
					op.local_Minvsigma2s[0].xdim,
					op.local_Minvsigma2s[0].ydim,
					baseMLO->mymodel.PPref[exp_iclass].data.xdim,
					baseMLO->mymodel.PPref[exp_iclass].data.ydim,
					baseMLO->mymodel.PPref[exp_iclass].data.zdim,
					baseMLO->mymodel.PPref[exp_iclass].data.yinit,
					baseMLO->mymodel.PPref[exp_iclass].data.zinit);
			model_real.free_device();
			model_imag.free_device();
			eulers.free();
			CUDA_CPU_TOC("generateModelProjections_wavg");
		}
		CUDA_CPU_TOC("projection_2");

		CudaGlobalPtr<FLOAT> wavgs_real(orientation_num * image_size);
		wavgs_real.device_alloc();
		wavgs_real.device_init(0);
		CudaGlobalPtr<FLOAT> wavgs_imag(orientation_num * image_size);
		wavgs_imag.device_alloc();
		wavgs_imag.device_init(0);
		CudaGlobalPtr<FLOAT> Fweights(orientation_num * image_size);
		Fweights.device_alloc();
		Fweights.device_init(0);

		/*=======================================================================================
										  PARTICLE ITERATION
		=======================================================================================*/

		/// Now that reference projection has been made loop over all particles inside this ori_particle
		for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
		{
			long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
			int group_id = baseMLO->mydata.getGroupId(part_id);

			double myprior_x, myprior_y, myprior_z;
			double old_offset_x = XX(op.old_offset[ipart]);
			double old_offset_y = YY(op.old_offset[ipart]);
			double old_offset_z;

			if (baseMLO->mymodel.ref_dim == 2)
			{
				myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
			}
			else
			{
				myprior_x = XX(op.prior[ipart]);
				myprior_y = YY(op.prior[ipart]);
				if (baseMLO->mymodel.data_dim == 3)
				{
					myprior_z = ZZ(op.prior[ipart]);
					old_offset_z = ZZ(op.old_offset[ipart]);
				}
			}


			/*======================================================
								 TRANSLATIONS
			======================================================*/

			CUDA_CPU_TIC("translation_2");

			CudaGlobalPtr<FLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans);
			CudaGlobalPtr<FLOAT> Fimgs_imag(Fimgs_real.size);
			CudaGlobalPtr<FLOAT> Fimgs_nomask_real(Fimgs_real.size);
			CudaGlobalPtr<FLOAT> Fimgs_nomask_imag(Fimgs_real.size);

			std::vector< long unsigned > iover_transes, itranses, ihiddens;

			long unsigned translation_num = imageTranslation(
					Fimgs_real,
					Fimgs_imag,
					Fimgs_nomask_real,
					Fimgs_nomask_imag,
					sp.itrans_min,
					sp.itrans_max,
					baseMLO->adaptive_oversampling ,
					baseMLO->sampling,
					oversampled_translations_x,
					oversampled_translations_y,
					oversampled_translations_z,
					sp.nr_oversampled_trans,
					baseMLO->global_fftshifts_ab_current,
					baseMLO->global_fftshifts_ab2_current,
					op.local_Fimgs_shifted[ipart],
					op.local_Fimgs_shifted_nomask[ipart],
					iover_transes,
					itranses,
					ihiddens,
					image_size);

			Fimgs_real.device_alloc();
			Fimgs_real.cp_to_device();
			Fimgs_imag.device_alloc();
			Fimgs_imag.cp_to_device();

			Fimgs_nomask_real.device_alloc();
			Fimgs_nomask_real.cp_to_device();
			Fimgs_nomask_imag.device_alloc();
			Fimgs_nomask_imag.cp_to_device();

			CUDA_CPU_TOC("translation_2");


			/*======================================================
					            	SCALE
			======================================================*/

			CUDA_CPU_TIC("scale_ctf");
			FLOAT part_scale(1.);

			if (baseMLO->do_scale_correction)
			{
				part_scale = baseMLO->mymodel.scale_correction[group_id];
				if (part_scale > 10000.)
				{
					std::cerr << " rlnMicrographScaleCorrection= " << part_scale << " group= " << group_id + 1 << std::endl;
					REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
				}
				else if (part_scale < 0.001)
				{
					if (!have_warned_small_scale)
					{
						std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << part_scale <<
								"); Use larger groups for more stable scale estimates." << std::endl;
						have_warned_small_scale = true;
					}
					part_scale = 0.001;
				}
			}

			CudaGlobalPtr<FLOAT> ctfs(image_size); //TODO Same size for all iparts, should be allocated once
			ctfs.device_alloc();

			if (baseMLO->do_ctf_correction)
			{
				for (unsigned i = 0; i < image_size; i++)
					ctfs[i] = (FLOAT) op.local_Fctfs[ipart].data[i] * part_scale;
			}
			else //TODO should be handled by memset
				for (unsigned i = 0; i < image_size; i++)
					ctfs[i] = part_scale;

			ctfs.cp_to_device();
			CUDA_CPU_TOC("scale_ctf");

			/*======================================================
					            MAP WEIGHTS
			======================================================*/

			CUDA_CPU_TIC("map");
			CudaGlobalPtr<FLOAT> sorted_weights(orientation_num * translation_num);

			mapWeights(
					sorted_weights,
					orientation_num,
					translation_num,
					baseMLO->sampling,
					ipart,
					iover_transes,
					ihiddens,
					iorientclasses,
					iover_rots,
					op.Mweight,
					sp.current_oversampling,
					sp.nr_trans);

			sorted_weights.device_alloc();
			sorted_weights.cp_to_device();
			sorted_weights.free_host();

			CUDA_CPU_TOC("map");

			/*======================================================
								KERNEL CALL
			======================================================*/

			// The below allocations are kept outside runWavgKernel(...) in case we decide to make them global.
			CudaGlobalPtr<FLOAT> Minvsigma2s(image_size); //TODO Same size for all iparts, should be allocated once
			Minvsigma2s.device_alloc();
			CudaGlobalPtr<FLOAT> wdiff2s_parts(orientation_num * image_size); //TODO Almost same size for all iparts, should be allocated once
			wdiff2s_parts.device_alloc();

			if (baseMLO->do_map)
				for (unsigned i = 0; i < image_size; i++)
					Minvsigma2s[i] = op.local_Minvsigma2s[ipart].data[i];
			else //TODO should be handled by memset
				for (unsigned i = 0; i < image_size; i++)
					Minvsigma2s[i] = 1;

			Minvsigma2s.cp_to_device();
			
#if !defined(CUDA_DOUBLE_PRECISION)
			if(do_combineProjAndWavg)
			{
				runProjAndWavgKernel(
						model_real,
						model_imag,
						eulers,
						Fimgs_real,
						Fimgs_imag,
						Fimgs_nomask_real,
						Fimgs_nomask_imag,
						sorted_weights,
						ctfs,
						Minvsigma2s,
						wdiff2s_parts,
						wavgs_real,
						wavgs_imag,
						Fweights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						image_size,
						XMIPP_MIN(baseMLO->mymodel.PPref[exp_iclass].r_max, op.local_Minvsigma2s[0].xdim - 1),
						ipart,
						group_id,
						exp_iclass
						);
			}
			else
#endif
			{
				runWavgKernel(
						Frefs_real,
						Frefs_imag,
						Fimgs_real,
						Fimgs_imag,
						Fimgs_nomask_real,
						Fimgs_nomask_imag,
						sorted_weights,
						ctfs,
						Minvsigma2s,
						wdiff2s_parts,
						wavgs_real,
						wavgs_imag,
						Fweights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass);

				Frefs_real.free_device();
				Frefs_imag.free_device();
			}
//			wdiff2s_parts.cp_to_host();
//			for (long unsigned k = 0; k < 100; k++)
//			{
//				std::cerr << wdiff2s_parts[k] << std::endl;
//			}
			/*======================================================
								COLLECT DATA
			======================================================*/

			CUDA_CPU_TIC("reduce_wdiff2s");
			// reduction_block_num = the highest possible power of two that covers more than or exactly half of all images to be reduced
			int num_reductions = (int)floor(log2((float)orientation_num));
			int reduction_block_num = pow(2,num_reductions);
			if(reduction_block_num==orientation_num) // (possibly) very special case where orientation_num is a power of 2
				reduction_block_num /= 2;

			CUDA_GPU_TIC("cuda_kernels_reduce_wdiff2s");
			unsigned orient1, orient2;
			for(int k=reduction_block_num; k>=1; k/=2) //invoke kernel repeatedly until all images have been stacked into the first image position
			{
				if(k>65535)
				{
					orient1 = ceil(sqrt(k));
					orient2 = orient1 + (orient1 % 2);  // For some reason the "optimal" values in the METADATA ooutput is sensitive to the choice of block-grid dims,
					orient1 +=          (orient1 % 2);  // and seems to work properly only when even numbers are used. // TODO examine why
				}
				else
				{
					orient1 = k;
					orient2 = 1;
				}
				dim3 grid_dim_wd(orient1,orient2);
				 // TODO **OF VERY LITTLE IMPORTANCE**  One block treating just 2 images is a very innefficient amount of loads per store
				cuda_kernel_reduce_wdiff2s<<<grid_dim_wd,BLOCK_SIZE>>>(~wdiff2s_parts,orientation_num,image_size,k);
			}
			CUDA_GPU_TOC("cuda_kernels_reduce_wdiff2s");

			wdiff2s_parts.size = image_size; //temporarily set the size to the single image we have now reduced, to not copy more than necessary
			wdiff2s_parts.cp_to_host();
			wdiff2s_parts.size = orientation_num * image_size;
			wdiff2s_parts.free_device();

			for (long int j = 0; j < image_size; j++)
			{
				int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, j);
				if (ires > -1)
				{
					thr_wsum_sigma2_noise[group_id].data[ires] += (double) wdiff2s_parts[j];
					exp_wsum_norm_correction[ipart] += (double) wdiff2s_parts[j];
				}
			}

			wdiff2s_parts.free_host();

			CUDA_CPU_TOC("reduce_wdiff2s");

			CUDA_CPU_TIC("collect_data_2");
			CUDA_CPU_TIC("collect_data_2_pre_kernel");
			//TODO should be replaced with loop over pairs of projections and translations (like in the getAllSquaredDifferences-function)

			std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;

			int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;
			CudaGlobalPtr<FLOAT >  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim]),
											sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples);
			int nr_transes = sp.nr_trans*sp.nr_oversampled_trans;
			CudaGlobalPtr<FLOAT>     oo_otrans_x(nr_transes); // old_offset_oversampled_trans_x
			CudaGlobalPtr<FLOAT>     oo_otrans_y(nr_transes);
			CudaGlobalPtr<FLOAT> myp_oo_otrans_x2y2z2(nr_transes); // my_prior_old_offs....x^2*y^2*z^2

			//Pregenerate oversampled translation objects for kernel-call
			for (long int itrans = 0, iitrans = 0; itrans < sp.nr_trans; itrans++)
			{
				baseMLO->sampling.getTranslations(itrans, baseMLO->adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
				{
					oo_otrans_x[iitrans] = old_offset_x + oversampled_translations_x[iover_trans];
					oo_otrans_y[iitrans] = old_offset_y + oversampled_translations_y[iover_trans];
					double diffx = myprior_x - oo_otrans_x[iitrans];
					double diffy = myprior_y - oo_otrans_y[iitrans];
					if (baseMLO->mymodel.data_dim == 3)
					{
						double diffz = myprior_z - (old_offset_z + oversampled_translations_z[iover_trans]);
						myp_oo_otrans_x2y2z2[iitrans] = diffx*diffx + diffy*diffy + diffz*diffz ;
					}
					else
					{
						myp_oo_otrans_x2y2z2[iitrans] = diffx*diffx + diffy*diffy ;
					}
				}
			}

			Mweight.device_alloc();
			Mweight.cp_to_device();
			oo_otrans_x.device_alloc();
			oo_otrans_x.cp_to_device();
			oo_otrans_y.device_alloc();
			oo_otrans_y.cp_to_device();
			myp_oo_otrans_x2y2z2.device_alloc();
			myp_oo_otrans_x2y2z2.cp_to_device();

			CudaGlobalPtr<FLOAT>                      p_weights(orientation_num);
			CudaGlobalPtr<FLOAT> p_thr_wsum_prior_offsetx_class(orientation_num);
			CudaGlobalPtr<FLOAT> p_thr_wsum_prior_offsety_class(orientation_num);
			CudaGlobalPtr<FLOAT>       p_thr_wsum_sigma2_offset(orientation_num);

			p_weights.device_alloc();
			p_thr_wsum_prior_offsetx_class.device_alloc();
			p_thr_wsum_prior_offsety_class.device_alloc();
			p_thr_wsum_sigma2_offset.device_alloc();

			dim3 grid_dim_collect2(sp.nr_dir, sp.nr_psi);
			CUDA_CPU_TOC("collect_data_2_pre_kernel");

			cuda_kernel_collect2<<<grid_dim_collect2,SUM_BLOCK_SIZE>>>(
					~oo_otrans_x,          // otrans-size -> make const
					~oo_otrans_y,          // otrans-size -> make const
					~myp_oo_otrans_x2y2z2, // otrans-size -> make const
					~Mweight,
					(FLOAT)op.significant_weight[ipart],
					(FLOAT)op.sum_weight[ipart],
					sp.nr_trans,
					sp.nr_oversampled_trans,
					sp.nr_oversampled_rot,
					oversamples,
					(baseMLO->do_skip_align || baseMLO->do_skip_rotate ),
					~p_weights,
					~p_thr_wsum_prior_offsetx_class,
					~p_thr_wsum_prior_offsety_class,
					~p_thr_wsum_sigma2_offset
				   );
			HANDLE_ERROR(cudaDeviceSynchronize());

			// TODO further reduce the below 4 arrays while data is still on gpu
			p_weights.cp_to_host();
			p_thr_wsum_prior_offsetx_class.cp_to_host();
			p_thr_wsum_prior_offsety_class.cp_to_host();
			p_thr_wsum_sigma2_offset.cp_to_host();

			thr_wsum_sigma2_offset = 0.0;
			int iorient = 0;
			for (long int idir = 0; idir < sp.nr_dir; idir++)
			{
				for (long int ipsi = 0; ipsi < sp.nr_psi; ipsi++, iorient++)
				{
					long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;
					// Only proceed if any of the particles had any significant coarsely sampled translation

					if (baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
					{
						long int mydir;
						if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
							mydir = idir;
						else
							mydir = op.pointer_dir_nonzeroprior[idir];

						// store partials according to indices of the relevant dimension
						DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += p_weights[iorient];
						thr_sumw_group[group_id]                 						+= p_weights[iorient];
						thr_wsum_pdf_class[exp_iclass]           						+= p_weights[iorient];
						thr_wsum_sigma2_offset                   						+= p_thr_wsum_sigma2_offset[iorient];

						if (baseMLO->mymodel.ref_dim == 2)
						{
							thr_wsum_prior_offsetx_class[exp_iclass] 	+= p_thr_wsum_prior_offsetx_class[iorient];
							thr_wsum_prior_offsety_class[exp_iclass] 	+= p_thr_wsum_prior_offsety_class[iorient];
						}
					}
				}
			}

			CUDA_CPU_TIC("collect_data_2_post_kernel");
			Mweight.free_device();
			p_weights.free();
			p_thr_wsum_sigma2_offset.free();
			p_thr_wsum_prior_offsetx_class.free();
			p_thr_wsum_prior_offsety_class.free();

			oo_otrans_y.free();
			oo_otrans_x.free();
			myp_oo_otrans_x2y2z2.free();

			//Get index of max element using GPU-tool thrust
			Indices max_index;
			max_index.fineIdx = thrust::max_element(&DIRECT_A2D_ELEM(op.Mweight, ipart, 0),&DIRECT_A2D_ELEM(op.Mweight, ipart+1, 0)) - &DIRECT_A2D_ELEM(op.Mweight, ipart, 0);
			op.max_weight[ipart] = DIRECT_A2D_ELEM(op.Mweight, ipart, max_index.fineIdx);
			max_index.fineIndexToFineIndices(sp); // set partial indices corresponding to the found max_index, to be used below

			baseMLO->sampling.getTranslations(max_index.itrans, baseMLO->adaptive_oversampling,
					oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
			baseMLO->sampling.getOrientations(max_index.idir, max_index.ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
					op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);
			double rot = oversampled_rot[max_index.ioverrot];
			double tilt = oversampled_tilt[max_index.ioverrot];
			double psi = oversampled_psi[max_index.ioverrot];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT) = rot;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT) = tilt;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PSI) = psi;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF) = XX(op.old_offset[ipart]) + oversampled_translations_x[max_index.iovertrans];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF) = YY(op.old_offset[ipart]) + oversampled_translations_y[max_index.iovertrans];
			if (baseMLO->mymodel.data_dim == 3)
				DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF) = ZZ(op.old_offset[ipart]) + oversampled_translations_z[max_index.iovertrans];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CLASS) = (double)max_index.iclass + 1;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX) = op.max_weight[ipart]/op.sum_weight[ipart];
			CUDA_CPU_TOC("collect_data_2_post_kernel");
			CUDA_CPU_TOC("collect_data_2");

		} // end loop ipart

		/*=======================================================================================
										   BACKPROJECTION
		=======================================================================================*/

		CUDA_CPU_TIC("backprojection");

		CudaGlobalPtr<FLOAT> bp_model_real(baseMLO->wsum_model.BPref[exp_iclass].data.nzyxdim);
		bp_model_real.device_alloc();
		bp_model_real.device_init(0);
		CudaGlobalPtr<FLOAT> bp_model_imag(bp_model_real.size);
		bp_model_imag.device_alloc();
		bp_model_imag.device_init(0);
		CudaGlobalPtr<FLOAT> bp_weight(bp_model_real.size);
		bp_weight.device_alloc();
		bp_weight.device_init(0);


		CudaGlobalPtr<FLOAT> bp_eulers(9 * orientation_num);

		FLOAT padding_factor = baseMLO->wsum_model.BPref[exp_iclass].padding_factor;

		generateEulerMatrices(
				1/padding_factor, //Why squared scale factor is given in backprojection
				rots,
				tilts,
				psis,
				bp_eulers,
				IS_NOT_INV);

		bp_eulers.device_alloc();
	    bp_eulers.cp_to_device();
	    bp_eulers.free_host();


	    runBackprojectKernel(
				wavgs_real,
				wavgs_imag,
				Fweights,
				bp_eulers,
				bp_model_real,
				bp_model_imag,
				bp_weight,
				baseMLO->wsum_model.BPref[exp_iclass].r_max,
				padding_factor * padding_factor,
				image_size,
				orientation_num,
				op.local_Minvsigma2s[0].xdim,
				op.local_Minvsigma2s[0].ydim,
				baseMLO->wsum_model.BPref[exp_iclass].data.xdim,
				baseMLO->wsum_model.BPref[exp_iclass].data.ydim,
				baseMLO->wsum_model.BPref[exp_iclass].data.zdim,
				baseMLO->wsum_model.BPref[exp_iclass].data.yinit,
				baseMLO->wsum_model.BPref[exp_iclass].data.zinit);

		bp_model_real.cp_to_host();
		bp_model_imag.cp_to_host();
		bp_weight.cp_to_host();

		HANDLE_ERROR(cudaDeviceSynchronize()); //TODO Optimize concurrency

		bp_model_real.free_device();
		bp_model_imag.free_device();
		bp_weight.free_device();

//#define PRINT_BACKPROJECTION_RESULTS
#ifdef PRINT_BACKPROJECTION_RESULTS

		FILE *fPtr1 = fopen("gpu_backproj_values.dat","w");
		for (unsigned i = 0; i < bp_model_real.size; i ++)
			fprintf(fPtr1, "%.1e %.1e\n", bp_model_real[i], bp_model_imag[i]);
		fclose(fPtr1);

		FILE *fPtr2 = fopen("gpu_backproj_weights.dat","w");
		for (unsigned i = 0; i < bp_weight.size; i ++)
			fprintf(fPtr2, "%.1e\n", bp_weight[i]);
		fclose(fPtr2);

		wavgs_real.cp_to_host();
		wavgs_imag.cp_to_host();
		Fweights.cp_to_host();

		for (long int i = 0; i < orientation_num; i++)
		{
			Euler_angles2matrix(rots[i], tilts[i], psis[i], A);

			for (unsigned j = 0; j < image_size; j++)
			{
				Fimg.data[j].real = (double) wavgs_real[i * image_size + j];
				Fimg.data[j].imag = (double) wavgs_imag[i * image_size + j];
				Fweight.data[j] = (double) Fweights[i * image_size + j];
			}

			int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
			pthread_mutex_lock(&global_mutex2[my_mutex]);
			(baseMLO->wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
			pthread_mutex_unlock(&global_mutex2[my_mutex]);

		}

		FILE *fPtr3 = fopen("cpu_backproj_values.dat","w");
		for (unsigned i = 0; i < (baseMLO->wsum_model.BPref[exp_iclass]).data.nzyxdim; i ++)
			fprintf(fPtr3, "%.1e %.1e\n", (baseMLO->wsum_model.BPref[exp_iclass]).data.data[i].real, (baseMLO->wsum_model.BPref[exp_iclass]).data.data[i].imag);
		fclose(fPtr3);

		FILE *fPtr4 = fopen("cpu_backproj_weights.dat","w");
		for (unsigned i = 0; i < (baseMLO->wsum_model.BPref[exp_iclass]).data.nzyxdim; i ++)
			fprintf(fPtr4, "%.1e\n", (baseMLO->wsum_model.BPref[exp_iclass]).bp_weight.data[i]);
		fclose(fPtr4);

		exit(0);
#endif

		Fweights.free();
		wavgs_real.free();
		wavgs_imag.free();

		int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
		pthread_mutex_lock(&global_mutex2[my_mutex]);

		for (long unsigned i = 0; i < bp_model_real.size; i++)
		{
			baseMLO->wsum_model.BPref[exp_iclass].data.data[i].real += bp_model_real[i];
			baseMLO->wsum_model.BPref[exp_iclass].data.data[i].imag += bp_model_imag[i];
			baseMLO->wsum_model.BPref[exp_iclass].weight.data[i] += bp_weight[i];
		}

		pthread_mutex_unlock(&global_mutex2[my_mutex]);

		CUDA_CPU_TOC("backprojection");

	} // end loop iclass

	CUDA_CPU_TIC("store_post_gpu");

	// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
	// Also calculate dLL for each particle and store in metadata
	// loop over all particles inside this ori_particle
	double thr_avg_norm_correction = 0.;
	double thr_sum_dLL = 0., thr_sum_Pmax = 0.;
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);

		// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
		for (int ires = baseMLO->mymodel.current_size/2 + 1; ires < baseMLO->mymodel.ori_size/2 + 1; ires++)
		{
			DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
			// Also extend the weighted sum of the norm_correction
			exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
		}

		// Store norm_correction
		// Multiply by old value because the old norm_correction term was already applied to the image
		if (baseMLO->do_norm_correction)
		{
			double old_norm_correction = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM);
			old_norm_correction /= baseMLO->mymodel.avg_norm_correction;
			// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
			// The variance of the total image (on which one normalizes) is twice this value!
			double normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
			thr_avg_norm_correction += normcorr;
			// Now set the new norm_correction in the relevant position of exp_metadata
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) = normcorr;


			// Print warning for strange norm-correction values
			if (!(baseMLO->iter == 1 && baseMLO->do_firstiter_cc) && DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) > 10.)
			{
				std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) << " for particle " << part_id << " in group " << group_id + 1 << "; Are your groups large enough?" << std::endl;
			}

		}

		// Store weighted sums for scale_correction
		if (baseMLO->do_scale_correction)
		{
			// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
			exp_wsum_scale_correction_XA[ipart] /= baseMLO->mymodel.scale_correction[group_id];
			exp_wsum_scale_correction_AA[ipart] /= baseMLO->mymodel.scale_correction[group_id] * baseMLO->mymodel.scale_correction[group_id];

			thr_wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
			thr_wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
		}

		// Calculate DLL for each particle
		double logsigma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine)
		{
			int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, n);
			// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
			// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
			if (ires > 0)
				logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], ires));
		}
		if (op.sum_weight[ipart]==0)
		{
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			std::cerr << " op.min_diff2[ipart]= " << op.min_diff2[ipart] << std::endl;
			std::cerr << " logsigma2= " << logsigma2 << std::endl;
			int group_id = baseMLO->mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.scale_correction[group_id]= " << baseMLO->mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_significant_weight[ipart]= " << op.significant_weight[ipart] << std::endl;
			std::cerr << " exp_max_weight[ipart]= " << op.max_weight[ipart] << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << baseMLO->mymodel.sigma2_noise[group_id] << std::endl;
			REPORT_ERROR("ERROR: op.sum_weight[ipart]==0");
		}
		double dLL;
		if ((baseMLO->iter==1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
			dLL = -op.min_diff2[ipart];
		else
			dLL = log(op.sum_weight[ipart]) - op.min_diff2[ipart] - logsigma2;

		// Store dLL of each image in the output array, and keep track of total sum
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_DLL) = dLL;
		thr_sum_dLL += dLL;

		// Also store sum of Pmax
		thr_sum_Pmax += DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX);

	}

	// Now, inside a global_mutex, update the other weighted sums among all threads
	if (!baseMLO->do_skip_maximization)
	{
		pthread_mutex_lock(&global_mutex);
		for (int n = 0; n < baseMLO->mymodel.nr_groups; n++)
		{
			baseMLO->wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			baseMLO->wsum_model.sumw_group[n] += thr_sumw_group[n];
			if (baseMLO->do_scale_correction)
			{
				baseMLO->wsum_model.wsum_signal_product_spectra[n] += thr_wsum_signal_product_spectra[n];
				baseMLO->wsum_model.wsum_reference_power_spectra[n] += thr_wsum_reference_power_spectra[n];
			}
		}
		for (int n = 0; n < baseMLO->mymodel.nr_classes; n++)
		{
			baseMLO->wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
			if (baseMLO->mymodel.ref_dim == 2)
			{
				XX(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
				YY(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
			}

			if (!(baseMLO->do_skip_align || baseMLO->do_skip_rotate) )
				baseMLO->wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		baseMLO->wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
		if (baseMLO->do_norm_correction)
			baseMLO->wsum_model.avg_norm_correction += thr_avg_norm_correction;
		baseMLO->wsum_model.LL += thr_sum_dLL;
		baseMLO->wsum_model.ave_Pmax += thr_sum_Pmax;
		pthread_mutex_unlock(&global_mutex);
	} // end if !do_skip_maximization

	CUDA_CPU_TOC("store_post_gpu");
}
