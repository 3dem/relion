#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <vector>
#include <iostream>
#include "src/acc/cuda/cuda_projector.h"
#include "src/acc/cuda/cuda_projector.cuh"
#include "src/acc/cuda/cuda_projector_plan.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_helper_functions.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"
#include "src/acc/acc_ptr.h"
#include "src/acc/cuda/cuda_ml_optimiser.h"
#include "src/complex.h"
#include "src/helix.h"
#include "src/error.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <signal.h>
#include <map>

#ifdef CUDA_FORCESTL
#include "src/acc/cuda/cuda_utils_stl.cuh"
#else
#include "src/acc/cuda/cuda_utils_cub.cuh"
#endif

#include "src/acc/acc_ml_optimiser_impl.h"


size_t MlDeviceBundle::checkFixedSizedObjects(int shares)
{
	int devCount;
	size_t BoxLimit;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	if(device_id >= devCount)
		CRITICAL(ERR_GPUID);

	HANDLE_ERROR(cudaSetDevice(device_id));

	size_t free(0), total(0);
	DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
	float margin(1.05);
	BoxLimit = pow(free/(margin*2.5*sizeof(XFLOAT)*((float)shares)),(1/3.0)) / (2.0);
	size_t BytesNeeded = ((float)shares)*margin*2.5*sizeof(XFLOAT)*pow((baseMLO->mymodel.ori_size*2),3);

	return(BoxLimit);
}
void MlDeviceBundle::setupFixedSizedObjects()
{
	unsigned nr_classes = baseMLO->mymodel.nr_classes;

	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(cudaSetDevice(device_id));

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR)
		generateProjectionPlanOnTheFly = true;
	else
		generateProjectionPlanOnTheFly = false;

	// clear() called on std::vector appears to set size=0, even if we have an explicit
	// destructor for each member, so we need to set the size to what is was before
	cudaProjectors.resize(nr_classes);
	cudaBackprojectors.resize(nr_classes);

	/*======================================================
	              PROJECTOR AND BACKPROJECTOR
	======================================================*/

	//Loop over classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		cudaProjectors[iclass].setMdlDim(
				baseMLO->mymodel.PPref[iclass].data.xdim,
				baseMLO->mymodel.PPref[iclass].data.ydim,
				baseMLO->mymodel.PPref[iclass].data.zdim,
				baseMLO->mymodel.PPref[iclass].data.yinit,
				baseMLO->mymodel.PPref[iclass].data.zinit,
				baseMLO->mymodel.PPref[iclass].r_max,
				baseMLO->mymodel.PPref[iclass].padding_factor);

		cudaProjectors[iclass].initMdl(baseMLO->mymodel.PPref[iclass].data.data);

		cudaBackprojectors[iclass].setMdlDim(
				baseMLO->wsum_model.BPref[iclass].data.xdim,
				baseMLO->wsum_model.BPref[iclass].data.ydim,
				baseMLO->wsum_model.BPref[iclass].data.zdim,
				baseMLO->wsum_model.BPref[iclass].data.yinit,
				baseMLO->wsum_model.BPref[iclass].data.zinit,
				baseMLO->wsum_model.BPref[iclass].r_max,
				baseMLO->wsum_model.BPref[iclass].padding_factor);

		cudaBackprojectors[iclass].initMdl();
	}

	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/

	int memAlignmentSize;
	cudaDeviceGetAttribute ( &memAlignmentSize, cudaDevAttrTextureAlignment, device_id );
	allocator = new CudaCustomAllocator(0, memAlignmentSize);
}

void MlDeviceBundle::setupTunableSizedObjects(size_t allocationSize)
{
	unsigned nr_classes = baseMLO->mymodel.nr_classes;
	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(cudaSetDevice(device_id));

	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/
#ifdef DEBUG_CUDA
	printf("DEBUG: Total GPU allocation size set to %zu MB on device id %d.\n", allocationSize / (1000*1000), device_id);
#endif
#ifndef CUDA_NO_CUSTOM_ALLOCATION
	allocator->resize(allocationSize);
#endif


	/*======================================================
	                    PROJECTION PLAN
	======================================================*/

	coarseProjectionPlans.resize(nr_classes, allocator);

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		//If doing predefined projector plan at all and is this class significant
		if (!generateProjectionPlanOnTheFly && baseMLO->mymodel.pdf_class[iclass] > 0.)
		{
			std::vector<int> exp_pointer_dir_nonzeroprior;
			std::vector<int> exp_pointer_psi_nonzeroprior;
			std::vector<RFLOAT> exp_directions_prior;
			std::vector<RFLOAT> exp_psi_prior;

			long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
			long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
			long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

			coarseProjectionPlans[iclass].setup(
					baseMLO->sampling,
					exp_directions_prior,
					exp_psi_prior,
					exp_pointer_dir_nonzeroprior,
					exp_pointer_psi_nonzeroprior,
					NULL, //Mcoarse_significant
					baseMLO->mymodel.pdf_class,
					baseMLO->mymodel.pdf_direction,
					nr_idir,
					nr_ipsi,
					0, //idir_min
					nr_idir - 1, //idir_max
					0, //ipsi_min
					nr_ipsi - 1, //ipsi_max
					0, //itrans_min
					itrans_max,
					0, //current_oversampling
					1, //nr_oversampled_rot
					iclass,
					true, //coarse
					!IS_NOT_INV,
					baseMLO->do_skip_align,
					baseMLO->do_skip_rotate,
					baseMLO->mymodel.orientational_prior_mode
					);
		}
	}
};

void MlOptimiserCuda::resetData()
{
	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(cudaSetDevice(device_id));

	unsigned nr_classes = baseMLO->mymodel.nr_classes;

	classStreams.resize(nr_classes, 0);
	for (int i = 0; i < nr_classes; i++)
		HANDLE_ERROR(cudaStreamCreate(&classStreams[i])); //HANDLE_ERROR(cudaStreamCreateWithFlags(&classStreams[i],cudaStreamNonBlocking));

	transformer1.clear();
	transformer2.clear();

	failsafe_attempts = 0;
};

void MlOptimiserCuda::doThreadExpectationSomeParticles(int thread_id)
{
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_THR);
#endif
//	CTOC(cudaMLO->timer,"interParticle");

	int devCount;
	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		DEBUG_HANDLE_ERROR(cudaSetDevice(device_id));
	//std::cerr << " calling on device " << device_id << std::endl;
	//put mweight allocation here
	size_t first_ipart = 0, last_ipart = 0;

	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		CTIC(timer,"oneTask");
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			CTIC(timer,"oneParticle");
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_A);
#endif
			unsigned my_ori_particle = baseMLO->exp_my_first_ori_particle + ipart;
			SamplingParameters sp;
			sp.nr_particles = baseMLO->mydata.ori_particles[my_ori_particle].particles_id.size();

			OptimisationParamters op(sp.nr_particles, my_ori_particle);

			// In the first iteration, multiple seeds will be generated
			// A single random class is selected for each pool of images, and one does not marginalise over the orientations
			// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
			// If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

			// Decide which classes to integrate over (for random class assignment in 1st iteration)
			sp.iclass_min = 0;
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
		    		// Make sure random division is always the same with the same seed
					long int idx = my_ori_particle - baseMLO->exp_my_first_ori_particle;
					if (idx >= baseMLO->exp_random_class_some_particles.size())
						REPORT_ERROR("BUG: expectationOneParticle idx>random_class_some_particles.size()");
					sp.iclass_min = sp.iclass_max = baseMLO->exp_random_class_some_particles[idx];
				}
			}

			// Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
			for (long int iori = baseMLO->exp_my_first_ori_particle; iori <= baseMLO->exp_my_last_ori_particle; iori++)
			{
				if (iori == my_ori_particle) break;
				op.metadata_offset += baseMLO->mydata.ori_particles[iori].particles_id.size();
			}
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_A);
#endif
			CTIC(timer,"getFourierTransformsAndCtfs");
			getFourierTransformsAndCtfs<MlOptimiserCuda,ACC_CUDA>(my_ori_particle, op, sp, baseMLO, this);
			CTOC(timer,"getFourierTransformsAndCtfs");

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
			//int nr_sampling_passes = (baseMLO->adaptive_oversampling > 0) ? 2 : 1;
			// But on the gpu the data-structures are different between passes, so we need to make a symbolic pass to set the weights up for storeWS
			int nr_sampling_passes = 2;

			/// -- This is a iframe-indexed vector, each entry of which is a dense data-array. These are replacements to using
			//    Mweight in the sparse (Fine-sampled) pass, coarse is unused but created empty input for convert ( FIXME )
			std::vector <IndexedDataArray> CoarsePassWeights(1, devBundle->allocator) ,FinePassWeights(sp.nr_particles, devBundle->allocator);
			// -- This is a iframe-indexed vector, each entry of which is a class-indexed vector of masks, one for each
			//    class in FinePassWeights
			std::vector < std::vector <IndexedDataArrayMask> > FinePassClassMasks(sp.nr_particles, std::vector <IndexedDataArrayMask>(baseMLO->mymodel.nr_classes, devBundle->allocator));
			// -- This is a iframe-indexed vector, each entry of which is parameters used in the projection-operations *after* the
			//    coarse pass, declared here to keep scope to storeWS
			std::vector < ProjectionParams > FineProjectionData(sp.nr_particles, baseMLO->mymodel.nr_classes);

			std::vector < cudaStager<unsigned long> > stagerD2(sp.nr_particles,devBundle->allocator), stagerSWS(sp.nr_particles,devBundle->allocator);

			for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
			{
				CTIC(timer,"weightPass");
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_B);
#endif
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
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_B);
#endif

				op.min_diff2.resize(sp.nr_particles, 0);
				op.avg_diff2.resize(sp.nr_particles, 0);

				if (ipass == 0)
				{
					unsigned long weightsPerPart(baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);

					op.Mweight.resizeNoCp(1,1,sp.nr_particles, weightsPerPart);

					CudaGlobalPtr<XFLOAT> Mweight(devBundle->allocator);
					Mweight.setSize(sp.nr_particles * weightsPerPart);
					Mweight.setHstPtr(op.Mweight.data);
					Mweight.device_alloc();
					deviceInitValue<XFLOAT>(Mweight, -999.);
					Mweight.streamSync();

					CTIC(timer,"getAllSquaredDifferencesCoarse");
					getAllSquaredDifferencesCoarse<MlOptimiserCuda,ACC_CUDA>(ipass, op, sp, baseMLO, this, Mweight);
					CTOC(timer,"getAllSquaredDifferencesCoarse");

					try
					{
						CTIC(timer,"convertAllSquaredDifferencesToWeightsCoarse");
						convertAllSquaredDifferencesToWeights<MlOptimiserCuda,ACC_CUDA,XFLOAT>(ipass, op, sp, baseMLO, this, CoarsePassWeights, FinePassClassMasks, Mweight);
						CTOC(timer,"convertAllSquaredDifferencesToWeightsCoarse");
					}
					catch (RelionError XE)
					{
						getAllSquaredDifferencesCoarse<MlOptimiserCuda,ACC_CUDA>(ipass, op, sp, baseMLO, this, Mweight);
#ifndef CUDA_DOUBLE_PRECISION
						try {
							convertAllSquaredDifferencesToWeights<MlOptimiserCuda,ACC_CUDA,double>(ipass, op, sp, baseMLO, this, CoarsePassWeights, FinePassClassMasks, Mweight);
						}
						catch (RelionError XE)
#endif
						{
							if (failsafe_attempts > 40)
								CRITICAL(ERRNUMFAILSAFE);

							//Rerun in fail-safe mode
							convertAllSquaredDifferencesToWeights<MlOptimiserCuda,ACC_CUDA,XFLOAT>(ipass, op, sp, baseMLO, this, CoarsePassWeights, FinePassClassMasks, Mweight, true);
							std::cerr << std::endl << "WARNING: Exception (" << XE.msg << ") handled by switching to fail-safe mode." << std::endl;
							failsafe_attempts ++;
						}
					}
				}
				else
				{
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_D);
#endif
//					// -- go through all classes and generate projectionsetups for all classes - to be used in getASDF and storeWS below --
//					// the reason to do this globally is subtle - we want the orientation_num of all classes to estimate a largest possible
//					// weight-array, which would be insanely much larger than necessary if we had to assume the worst.
					for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
					{
						FineProjectionData[iframe].orientationNumAllClasses = 0;
						for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
						{
							if(exp_iclass>0)
								FineProjectionData[iframe].class_idx[exp_iclass] = FineProjectionData[iframe].rots.size();
							FineProjectionData[iframe].class_entries[exp_iclass] = 0;

							CTIC(timer,"generateProjectionSetup");
							FineProjectionData[iframe].orientationNumAllClasses += generateProjectionSetupFine(
									op,
									sp,
									baseMLO,
									exp_iclass,
									FineProjectionData[iframe]);
							CTOC(timer,"generateProjectionSetup");

						}
						//set a maximum possible size for all weights (to be reduced by significance-checks)
						FinePassWeights[iframe].setDataSize(FineProjectionData[iframe].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans);
						FinePassWeights[iframe].dual_alloc_all();
						stagerD2[iframe].size= 2*(FineProjectionData[iframe].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans);
						stagerD2[iframe].prepare();
					}
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_D);
#endif
//					printf("Allocator used space before 'getAllSquaredDifferencesFine': %.2f MiB\n", (float)devBundle->allocator->getTotalUsedSpace()/(1024.*1024.));

					CTIC(timer,"getAllSquaredDifferencesFine");
					getAllSquaredDifferencesFine<MlOptimiserCuda,ACC_CUDA>(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks, FineProjectionData, stagerD2);
					CTOC(timer,"getAllSquaredDifferencesFine");
					FinePassWeights[0].weights.cp_to_host();
					CudaGlobalPtr<XFLOAT> Mweight(devBundle->allocator); //DUMMY

					CTIC(timer,"convertAllSquaredDifferencesToWeightsFine");
					convertAllSquaredDifferencesToWeights<MlOptimiserCuda,ACC_CUDA,XFLOAT>(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks, Mweight);
					CTOC(timer,"convertAllSquaredDifferencesToWeightsFine");

				}

				CTOC(timer,"weightPass");
			}
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_E);
#endif

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;
			for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
			{
				stagerSWS[iframe].size= 2*(FineProjectionData[iframe].orientationNumAllClasses);
				stagerSWS[iframe].prepare();
			}
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_E);
#endif
			CTIC(timer,"storeWeightedSums");
			storeWeightedSums<MlOptimiserCuda,ACC_CUDA>(op, sp, baseMLO, this, FinePassWeights, FineProjectionData, FinePassClassMasks, stagerSWS);
			CTOC(timer,"storeWeightedSums");

			CTOC(timer,"oneParticle");
		}
		CTOC(timer,"oneTask");
	}

//	CTIC(cudaMLO->timer,"interParticle");
//	exit(0);

#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_THR);
#endif
}

