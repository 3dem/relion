#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/cuda/cuda_projector_plan.h"
#include "src/acc/cuda/cuda_projector.h"
#include "src/acc/cuda/cuda_backprojector.h"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#include <stack>
//#include <cufft.h>


#include "src/acc/acc_ml_optimiser.h"

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


class IndexedDataArrayMask
{
public:
	// indexes of job partition
	//   every element in jobOrigin    is a reference to point to a position in a IndexedDataArray.weights array where that job starts RELATIVE to firstPos
	//   every element in jobExtent    specifies the number of weights for that job
	CudaGlobalPtr<size_t> jobOrigin, jobExtent;

	size_t firstPos, lastPos; // positions in indexedDataArray data and index arrays to slice out
	size_t weightNum, jobNum; // number of weights and jobs this class

	inline
	 IndexedDataArrayMask(CudaCustomAllocator *allocator):
	jobOrigin(allocator),
	jobExtent(allocator),
	firstPos(),
	lastPos(),
	weightNum(),
	jobNum()
	{};

public:

	void setNumberOfJobs(size_t newSize)
	{
		jobNum=newSize;
		jobOrigin.setSize(newSize);
		jobExtent.setSize(newSize);
	}

	void setNumberOfWeights(size_t newSize)
	{
		weightNum=newSize;
	}

	inline
	 ~IndexedDataArrayMask()
	{
//		jobOrigin.free_host();
//		jobExtent.free_host();
	};
};

class IndexedDataArray
{
public:
	//actual data
	CudaGlobalPtr<XFLOAT> weights;

	// indexes with same length as data
	// -- basic indices ---------------------------------
	//     rot_id  = id of rot     = which of all POSSIBLE orientations                               this weight signifies
	//     rot_idx = index of rot  = which in the sequence of the determined significant orientations this weight signifies
	//   trans_id  = id of trans   = which of all POSSIBLE translations                               this weight signifies
	// -- special indices ---------------------------------
	//   ihidden_overs  =  mapping to MWeight-based indexing for compatibility
	CudaGlobalPtr<size_t> rot_id, rot_idx, trans_idx, ihidden_overs;

	inline
	 IndexedDataArray(CudaCustomAllocator *allocator):
		weights(allocator),
		rot_id(allocator),
		rot_idx(allocator),
		trans_idx(allocator),
		ihidden_overs(allocator)
	{};

	// constructor which takes a parent IndexedDataArray and a mask to create a child
	inline
	 IndexedDataArray(IndexedDataArray &parent, IndexedDataArrayMask &mask, CudaCustomAllocator *allocator):
		weights(		&(parent.weights.h_ptr[mask.firstPos])		,&(parent.weights.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		rot_id(			&(parent.rot_id.h_ptr[mask.firstPos])		,&(parent.rot_id.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		rot_idx(		&(parent.rot_idx.h_ptr[mask.firstPos])		,&(parent.rot_idx.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		trans_idx(		&(parent.trans_idx.h_ptr[mask.firstPos])	,&(parent.trans_idx.d_ptr[mask.firstPos])		,mask.weightNum, allocator),
		ihidden_overs(	&(parent.ihidden_overs.h_ptr[mask.firstPos]),&(parent.ihidden_overs.d_ptr[mask.firstPos])	,mask.weightNum, allocator)
	{
		weights.d_do_free=false;
		rot_id.d_do_free=false;
		rot_idx.d_do_free=false;
		trans_idx.d_do_free=false;
		ihidden_overs.d_do_free=false;

		weights.h_do_free=false;
		rot_id.h_do_free=false;
		rot_idx.h_do_free=false;
		trans_idx.h_do_free=false;
		ihidden_overs.h_do_free=false;
	};

public:

	void setDataSize(size_t newSize)
	{
		weights.setSize(newSize);
		rot_id.setSize(newSize);
		rot_idx.setSize(newSize);
		trans_idx.setSize(newSize);
		ihidden_overs.setSize(newSize);
	}

	void resize_host_all(size_t newSize)
	{
		weights.resize_host(newSize);
		rot_id.resize_host(newSize);
		rot_idx.resize_host(newSize);
		trans_idx.resize_host(newSize);
		ihidden_overs.resize_host(newSize);
	}

	void host_alloc_all()
	{
		weights.host_alloc();
		rot_id.host_alloc();
		rot_idx.host_alloc();
		trans_idx.host_alloc();
		ihidden_overs.host_alloc();
	}

	void device_alloc_all()
	{
		weights.device_alloc();
		rot_id.device_alloc();
		rot_idx.device_alloc();
		trans_idx.device_alloc();
		ihidden_overs.device_alloc();
	}

	void dual_alloc_all()
	{
		host_alloc_all();
		device_alloc_all();
	}
};


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
	std::vector< CudaProjectorPlan > coarseProjectionPlans;

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
