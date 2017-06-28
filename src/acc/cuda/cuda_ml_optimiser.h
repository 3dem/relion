#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cuda/cuda_projector.h"
#include "src/acc/cuda/cuda_backprojector.h"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#include <stack>
//#include <cufft.h>


#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ptr.h"

//#ifdef DEBUG_CUDA
//#define HANDLE_CUFFT_ERROR( err ) (CufftHandleError( err, __FILE__, __LINE__ ))
//#else
//#define HANDLE_CUFFT_ERROR( err ) (err) //Do nothing
//#endif
//static void CufftHandleError( cufftResult err, const char *file, int line )
//{
//    if (err != CUFFT_SUCCESS)
//    {
//        fprintf(stderr, "Cufft error in file '%s' in line %i : %s.\n",
//                __FILE__, __LINE__, "error" );
//		raise(SIGSEGV);
//    }
//}

/*
 * Bundle of device-objects
 */
class MlDeviceBundle
{
public:

	//The CUDA accelerated projector set
	std::vector< CudaProjector > cudaProjectors;

	//The CUDA accelerated back-projector set
	std::vector< CudaBackprojector > cudaBackprojectors;

	//Used for precalculations of projection setup
	CudaCustomAllocator *allocator;

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
		DEBUG_HANDLE_ERROR(cudaDeviceSynchronize());
	}


	~MlDeviceBundle()
	{
		cudaProjectors.clear();
		cudaBackprojectors.clear();
		coarseProjectionPlans.clear();
		//Delete this lastly
		delete allocator;
		HANDLE_ERROR(cudaSetDevice(device_id));
		HANDLE_ERROR(cudaDeviceReset());
	}

};

class MlOptimiserCuda
{
public:
	// transformer as holder for reuse of fftw_plans
	FourierTransformer transformer;

   //Class streams ( for concurrent scheduling of class-specific kernels)
	std::vector< cudaStream_t > classStreams;
	cudaError_t errorStatus;

	CudaFFT transformer1;
	CudaFFT transformer2;

	MlOptimiser *baseMLO;

	bool refIs3D;
	bool dataIs3D;

	int device_id;

	unsigned failsafe_attempts;

	MlDeviceBundle *devBundle;

#ifdef TIMING_FILES
	relion_timer timer;
#endif

	MlOptimiserCuda(MlOptimiser *baseMLOptimiser, MlDeviceBundle* bundle, const char * timing_fnm) :
			baseMLO(baseMLOptimiser),
			transformer1(cudaStreamPerThread, bundle->allocator, baseMLOptimiser->mymodel.data_dim),
			transformer2(cudaStreamPerThread, bundle->allocator, baseMLOptimiser->mymodel.data_dim),
			refIs3D(baseMLO->mymodel.ref_dim == 3),
			dataIs3D(baseMLO->mymodel.data_dim == 3),
			devBundle(bundle),
			device_id(bundle->device_id),
#ifdef TIMING_FILES
			timer(timing_fnm),
#endif
			errorStatus((cudaError_t)0),
			failsafe_attempts(0)
	{};

	void resetData();

	void doThreadExpectationSomeParticles(int thread_id);

	~MlOptimiserCuda()
	{
		for (int i = 0; i < classStreams.size(); i++)
			if (classStreams[i] != NULL)
				HANDLE_ERROR(cudaStreamDestroy(classStreams[i]));
	}

};

#endif
