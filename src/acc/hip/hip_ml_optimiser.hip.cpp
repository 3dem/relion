/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#undef ALTCPU
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <vector>
#include <iostream>
#include "src/ml_optimiser.h"
#include <hip/hip_runtime.h>
#include <hiprand/hiprand.h>
#include <hiprand/hiprand_kernel.h>

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/hip/hip_benchmark_utils.h"
#include "src/acc/hip/hip_kernels/helper.h"
#include "src/acc/hip/hip_kernels/diff2.h"
#include "src/acc/hip/hip_kernels/wavg.h"
#include "src/acc/hip/hip_mem_utils.h"
#include "src/acc/hip/hip_fft.h"
#include "src/acc/data_types.h"
#include "src/complex.h"
#include "src/helix.h"
#include "src/error.h"
#include <fstream>
#include "src/parallel.h"
#include <signal.h>
#include <map>

#ifdef HIP_FORCESTL
#include "src/acc/hip/hip_utils_stl.h"
#else
#include "src/acc/hip/hip_utils_cub.h"
#endif

#include "src/acc/utilities.h"
#include "src/acc/utilities_impl.h"

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/hip/hip_ml_optimiser.h"
#include "src/acc/acc_helper_functions.h"
#include "src/acc/acc_ml_optimiser_impl.h"

// -------------------------------  Some explicit template instantiations
template __global__ void HipKernels::hip_kernel_translate2D<XFLOAT>(XFLOAT *,
    XFLOAT*, int, int, int, int, int);

template __global__ void HipKernels::hip_kernel_translate3D<XFLOAT>(XFLOAT *,
    XFLOAT *, int, int, int, int, int, int, int);

template __global__ void hip_kernel_multi<XFLOAT>( XFLOAT *,
	XFLOAT *, XFLOAT, int);

template __global__ void HipKernels::hip_kernel_multi<XFLOAT>( XFLOAT *,
	XFLOAT, int);

template __global__ void hip_kernel_multi<XFLOAT>( XFLOAT *, XFLOAT *,
	XFLOAT *, XFLOAT, int);

// ----------------------------------------------------------------------

// High-level HIP objects

size_t MlDeviceBundle::checkFixedSizedObjects(int shares)
{
	int devCount;
	size_t BoxLimit;
	HANDLE_ERROR(hipGetDeviceCount(&devCount));
	if(device_id >= devCount)
		CRITICAL(ERR_GPUID);

	HANDLE_ERROR(hipSetDevice(device_id));

	size_t free(0), total(0);
	DEBUG_HANDLE_ERROR(hipMemGetInfo( &free, &total ));
	float margin(1.05);
	BoxLimit = pow(free/(margin*2.5*sizeof(XFLOAT)*((float)shares)),(1/3.0)) / ((float) baseMLO->mymodel.padding_factor);
	//size_t BytesNeeded = ((float)shares)*margin*2.5*sizeof(XFLOAT)*pow((baseMLO->mymodel.ori_size*baseMLO->mymodel.padding_factor),3);

	return(BoxLimit);
}
void MlDeviceBundle::setupFixedSizedObjects()
{
	int devCount;
	HANDLE_ERROR(hipGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(hipSetDevice(device_id));

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR || baseMLO->mydata.is_tomo)
		generateProjectionPlanOnTheFly = true;
	else
		generateProjectionPlanOnTheFly = false;

	unsigned nr_proj = baseMLO->mymodel.PPref.size();
	unsigned nr_bproj = baseMLO->wsum_model.BPref.size();

	projectors.resize(nr_proj);
	backprojectors.resize(nr_bproj);

	/*======================================================
	              PROJECTOR AND BACKPROJECTOR
	======================================================*/

	for (int imodel = 0; imodel < nr_proj; imodel++)
	{
		projectors[imodel].setMdlDim(
				baseMLO->mymodel.PPref[imodel].data.xdim,
				baseMLO->mymodel.PPref[imodel].data.ydim,
				baseMLO->mymodel.PPref[imodel].data.zdim,
				baseMLO->mymodel.PPref[imodel].data.yinit,
				baseMLO->mymodel.PPref[imodel].data.zinit,
				baseMLO->mymodel.PPref[imodel].r_max,
				baseMLO->mymodel.PPref[imodel].padding_factor);

		projectors[imodel].initMdl(baseMLO->mymodel.PPref[imodel].data.data);

	}

	for (int imodel = 0; imodel < nr_bproj; imodel++)
	{
		backprojectors[imodel].setMdlDim(
				baseMLO->wsum_model.BPref[imodel].data.xdim,
				baseMLO->wsum_model.BPref[imodel].data.ydim,
				baseMLO->wsum_model.BPref[imodel].data.zdim,
				baseMLO->wsum_model.BPref[imodel].data.yinit,
				baseMLO->wsum_model.BPref[imodel].data.zinit,
				baseMLO->wsum_model.BPref[imodel].r_max,
				baseMLO->wsum_model.BPref[imodel].padding_factor);

		backprojectors[imodel].initMdl();
	}

	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/

	int memAlignmentSize;
	hipDeviceGetAttribute ( &memAlignmentSize, hipDeviceAttributeTextureAlignment, device_id );
	allocator = new HipCustomAllocator(0, memAlignmentSize);
}

void MlDeviceBundle::setupTunableSizedObjects(size_t allocationSize)
{
	unsigned nr_models = baseMLO->mymodel.nr_classes;
	int devCount;
	HANDLE_ERROR(hipGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(hipSetDevice(device_id));

	/*======================================================
	                    CUSTOM ALLOCATOR
	======================================================*/
#ifdef DEBUG_HIP
	printf("DEBUG: Total GPU allocation size set to %zu MB on device id %d.\n", allocationSize / (1000*1000), device_id);
#endif
#ifndef HIP_NO_CUSTOM_ALLOCATION
	allocator->resize(allocationSize);
#endif


	/*======================================================
	                    PROJECTION PLAN
	======================================================*/

	coarseProjectionPlans.resize(nr_models, allocator);

	for (int iclass = 0; iclass < nr_models; iclass++)
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

void MlOptimiserHip::resetData()
{
	int devCount;
	HANDLE_ERROR(hipGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		HANDLE_ERROR(hipSetDevice(device_id));

	unsigned nr_classes = baseMLO->mymodel.nr_classes;

	classStreams.resize(nr_classes, 0);
	for (int i = 0; i < nr_classes; i++)
		HANDLE_ERROR(hipStreamCreate(&classStreams[i])); //HANDLE_ERROR(hipStreamCreateWithFlags(&classStreams[i],hipStreamNonBlocking));

	transformer1.clear();
	transformer2.clear();
};

void MlOptimiserHip::doThreadExpectationSomeParticles(int thread_id)
{
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_THR);
#endif
//	CTOC(hipMLO->timer,"interParticle");

	int devCount;
	HANDLE_ERROR(hipGetDeviceCount(&devCount));
	if(device_id >= devCount)
	{
		//std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		CRITICAL(ERR_GPUID);
	}
	else
		DEBUG_HANDLE_ERROR(hipSetDevice(device_id));
	//std::cerr << " calling on device " << device_id << std::endl;
	//put mweight allocation here
	size_t first_ipart = 0, last_ipart = 0;

	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		CTIC(timer,"oneTask");
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_A);
#endif

			AccPtrFactory ptrFactory(allocator, hipStreamPerThread);
            accDoExpectationOneParticle<MlOptimiserHip>(this, baseMLO->exp_my_first_part_id + ipart, thread_id, ptrFactory);

		}
		CTOC(timer,"oneTask");
	}

//	CTIC(hipMLO->timer,"interParticle");
//	exit(0);

#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_THR);
#endif
}
