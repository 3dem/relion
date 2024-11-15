/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef HIP_ML_OPTIMISER_H_
#define HIP_ML_OPTIMISER_H_

#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/hip/hip_mem_utils.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/hip/hip_fft.h"
#include "src/acc/hip/hip_benchmark_utils.h"
#include <stack>
//#include <hipfft.h>

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ptr.h"

class MlDeviceBundle
{
public:

	//The HIP accelerated projector set
	std::vector< AccProjector > projectors;

	//The HIP accelerated back-projector set
	std::vector< AccBackprojector > backprojectors;

	//Used for precalculations of projection setup
	HipCustomAllocator *allocator;

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;
	std::vector< AccProjectorPlan > coarseProjectionPlans;

	MlOptimiser *baseMLO;

	int device_id;

	int rank_shared_count;

	bool haveWarnedRefinementMem;

	MlDeviceBundle(MlOptimiser *baseMLOptimiser):
			baseMLO(baseMLOptimiser),
			generateProjectionPlanOnTheFly(false),
			rank_shared_count(1),
			device_id(-1),
			haveWarnedRefinementMem(false),
			allocator(NULL)
	{};

	void setDevice(int did)
	{
		device_id = did;
	}

	size_t checkFixedSizedObjects(int shares);
	void setupFixedSizedObjects();
	void setupTunableSizedObjects(size_t allocationSize);

	void syncAllBackprojects()
	{
		DEBUG_HANDLE_ERROR(hipDeviceSynchronize());
	}


	~MlDeviceBundle()
	{
		projectors.clear();
		backprojectors.clear();
		coarseProjectionPlans.clear();
		//Delete this lastly
		delete allocator;
	}

};
class MlOptimiserHip
{
public:
	// transformer as holder for reuse of fftw_plans
	FourierTransformer transformer;

   //Class streams ( for concurrent scheduling of class-specific kernels)
	std::vector< hipStream_t > classStreams;
	hipStream_t defaultStream = 0;
	hipError_t errorStatus;

	HipFFT transformer1;
	HipFFT transformer2;

	MlOptimiser *baseMLO;

	bool refIs3D;
	bool dataIs3D;
        bool shiftsIs3D;

	int device_id;

	MlDeviceBundle *bundle;

	//Used for precalculations of projection setup
	HipCustomAllocator *allocator;

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;


#ifdef TIMING_FILES
	relion_timer timer;
#endif

	MlOptimiserHip(MlOptimiser *baseMLOptimiser, MlDeviceBundle* bundle, const char * timing_fnm) :
			baseMLO(baseMLOptimiser),
			transformer1(hipStreamPerThread, bundle->allocator, baseMLOptimiser->mymodel.data_dim),
			transformer2(hipStreamPerThread, bundle->allocator, baseMLOptimiser->mymodel.data_dim),
			refIs3D(baseMLO->mymodel.ref_dim == 3),
			dataIs3D(baseMLO->mymodel.data_dim == 3),
			shiftsIs3D(baseMLO->mymodel.data_dim == 3 || baseMLO->mydata.is_tomo),
            bundle(bundle),
			device_id(bundle->device_id),
#ifdef TIMING_FILES
			timer(timing_fnm),
#endif
			errorStatus((hipError_t)0),
			allocator(bundle->allocator),
			generateProjectionPlanOnTheFly(bundle->generateProjectionPlanOnTheFly)
	{};

	void resetData();

	void doThreadExpectationSomeParticles(int thread_id);

	~MlOptimiserHip()
	{
		for (int i = 0; i < classStreams.size(); i++)
			if (classStreams[i] != NULL)
				HANDLE_ERROR(hipStreamDestroy(classStreams[i]));
	}

	HipCustomAllocator *getAllocator()
	{
		return (bundle->allocator);
	};

};

#endif
