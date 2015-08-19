#ifndef CUDA_UTILS_CUH_
#define CUDA_UTILS_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <stdio.h>
#include <signal.h>
#include <vector>
#include <src/gpu_utils/cub/cub.cuh>

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
__device__ inline XFLOAT cuda_atomic_add(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return __longlong_as_double(old);
}
#else
#define XFLOAT float
__device__ inline void cuda_atomic_add(float* address, float value)
{
  atomicAdd(address,value);
}
#endif

template <typename T>
inline static cub::KeyValuePair<int, T> getArgMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
	CudaGlobalPtr<cub::KeyValuePair<int, T> >  max_pair(1, ptr.getAllocator());
	max_pair.device_alloc();
	size_t temp_storage_size = 0;

	cub::DeviceReduce::ArgMax( NULL, temp_storage_size, ~ptr, ~max_pair, ptr.size);

	CudaCustomAllocator::Alloc* alloc = ptr.getAllocator()->alloc(temp_storage_size);

	cub::DeviceReduce::ArgMax( alloc->getPtr(), temp_storage_size, ~ptr, ~max_pair, ptr.size);

	max_pair.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(ptr.getStream()));

	ptr.getAllocator()->free(alloc);

	return max_pair[0];
}


class IndexedDataArrayMask
{
public:
	// indexes of job partition
	//   every element in jobOrigin    is a reference to point to a position in a IndexedDataArray.weights array where that job starts RELATIVE to firstPos
	//   every element in jobExtent    specifies the number of weights for that job
	CudaGlobalPtr<long unsigned> jobOrigin, jobExtent;

	long unsigned firstPos, lastPos; // positions in indexedDataArray data and index arrays to slice out
	long unsigned weightNum, jobNum; // number of weights and jobs this class

	inline
	__host__  IndexedDataArrayMask(CudaCustomAllocator *allocator):
	jobOrigin(allocator),
	jobExtent(allocator),
	firstPos(),
	lastPos(),
	weightNum(),
	jobNum()
	{};

public:

	__host__ void setNumberOfJobs(long int newSize)
	{
		jobNum=newSize;
		jobOrigin.size=newSize;
		jobExtent.size=newSize;
	}

	__host__ void setNumberOfWeights(long int newSize)
	{
		weightNum=newSize;
	}

	inline
	__host__  ~IndexedDataArrayMask()
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
	//   class_id  = id of class   = which of all POSSIBLE classes                                    this weight signifies
	// -- special indices ---------------------------------
	//   ihidden_overs  =  mapping to MWeight-based indexing for compatibility
	CudaGlobalPtr<long unsigned> rot_id, rot_idx, trans_idx, ihidden_overs, class_id;

	inline
	__host__  IndexedDataArray(CudaCustomAllocator *allocator):
		weights(allocator),
		rot_id(allocator),
		rot_idx(allocator),
		trans_idx(allocator),
		ihidden_overs(allocator),
		class_id(allocator)
	{};

	// constructor which takes a parent IndexedDataArray and a mask to create a child
	inline
	__host__  IndexedDataArray(IndexedDataArray &parent, IndexedDataArrayMask &mask, CudaCustomAllocator *allocator):
		weights(		&(parent.weights.h_ptr[mask.firstPos])		,&(parent.weights.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		rot_id(			&(parent.rot_id.h_ptr[mask.firstPos])		,&(parent.rot_id.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		rot_idx(		&(parent.rot_idx.h_ptr[mask.firstPos])		,&(parent.rot_idx.d_ptr[mask.firstPos])			,mask.weightNum, allocator),
		trans_idx(		&(parent.trans_idx.h_ptr[mask.firstPos])	,&(parent.trans_idx.d_ptr[mask.firstPos])		,mask.weightNum, allocator),
		ihidden_overs(	&(parent.ihidden_overs.h_ptr[mask.firstPos]),&(parent.ihidden_overs.d_ptr[mask.firstPos])	,mask.weightNum, allocator),
		class_id(		&(parent.class_id.h_ptr[mask.firstPos])		,&(parent.class_id.d_ptr[mask.firstPos])		,mask.weightNum, allocator)
	{
		weights.d_do_free=false;
		rot_id.d_do_free=false;
		rot_idx.d_do_free=false;
		trans_idx.d_do_free=false;
		ihidden_overs.d_do_free=false;
		class_id.d_do_free=false;

		weights.h_do_free=false;
		rot_id.h_do_free=false;
		rot_idx.h_do_free=false;
		trans_idx.h_do_free=false;
		ihidden_overs.h_do_free=false;
		class_id.h_do_free=false;
	};

public:

	__host__ void setDataSize(long int newSize)
	{
		weights.size=newSize;
		rot_id.size=newSize;
		rot_idx.size=newSize;
		trans_idx.size=newSize;
		ihidden_overs.size=newSize;
		class_id.size=newSize;
	}

	__host__ void dual_alloc_all()
	{
		weights.host_alloc();
		rot_id.host_alloc();
		rot_idx.host_alloc();
		trans_idx.host_alloc();
		ihidden_overs.host_alloc();
		class_id.host_alloc();
		//-----------------------
		weights.device_alloc();
		rot_id.device_alloc();
		rot_idx.device_alloc();
		trans_idx.device_alloc();
		ihidden_overs.device_alloc();
		class_id.device_alloc();
	}
};


class ProjectionParams
{

public:
	std::vector< long unsigned > orientation_num; 					// the number of significant orientation for each class
	long unsigned orientationNumAllClasses;							// sum of the above
	std::vector< double > rots, tilts, psis;
	std::vector< long unsigned > iorientclasses, iover_rots;

	// These are arrays which detial the number of entries in each class, and where each class starts.
	// NOTE: There is no information about which class each class_idx refers to, there is only
	// a distinction between different classes.
	std::vector< long unsigned > class_entries, class_idx;
	inline
	__host__ ProjectionParams():

		rots(),
		tilts(),
		psis(),
		iorientclasses(),
		iover_rots(),

		class_entries(),
		class_idx(),
		orientation_num(),
		orientationNumAllClasses(0)

	{};

	inline
	__host__ ProjectionParams(unsigned long classes):

		rots(),
		tilts(),
		psis(),
		iorientclasses(),
		iover_rots(),

		class_entries(classes),
		class_idx(classes),
		orientation_num(classes),
		orientationNumAllClasses(0)
	{
		class_idx[0]=0;
		class_entries[0]=0;
	};


	// constructor that slices out a part of a parent ProjectionParams, assumed to contain a single (partial or entire) class
	inline
	__host__ ProjectionParams(ProjectionParams &parent, unsigned long start, unsigned long end):
		rots(			&parent.rots[start],  			&parent.rots[end]),
		tilts(			&parent.tilts[start], 			&parent.tilts[end]),
		psis(			&parent.psis[start],  			&parent.psis[end]),
		iorientclasses( &parent.iorientclasses[start],  &parent.iorientclasses[end]),
		iover_rots(		&parent.iover_rots[start],  	&parent.iover_rots[end]),
		orientation_num(1),
		class_entries(1,end-start),
		class_idx(1,0) // NOTE: this is NOT the class, but rather where in these partial PrjParams to start, which is @ 0.
	{};

public:
	// Appends new values into the projection parameters for later use.
	// class_idx is used as such:
	// the n:th class (beginning with 0:th)
	// begins @ element class_idx[n]
	// ends   @ element class_idx[n]+class_entries[n]

	__host__ void pushBackAll(long unsigned iclass, double NEWrot,double NEWtilt ,double NEWpsi, long unsigned NEWiorientclasses,long unsigned NEWiover_rots)
	{
		// incremement the counter for this class
		class_entries[iclass]++;
		// and push a new entry
		rots.push_back(NEWrot);
		tilts.push_back(NEWtilt);
		psis.push_back(NEWpsi);
		iorientclasses.push_back(NEWiorientclasses);
		iover_rots.push_back(NEWiover_rots);
	}
};
#endif

