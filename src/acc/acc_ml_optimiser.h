#ifndef ACC_ML_OPTIMISER_H_
#define ACC_ML_OPTIMISER_H_
#if 0
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
#endif

#include "src/acc/acc_ptr.h"

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif


class SamplingParameters
{
public:
	unsigned long nr_dir,
	nr_psi,
	nr_trans,
	nr_oversampled_rot,
	nr_oversampled_trans,
	nr_particles,
	current_oversampling,
	current_image_size,
	iclass_min, iclass_max,
	idir_min, idir_max,
	ipsi_min, ipsi_max,
	itrans_min, itrans_max;
	std::string current_img;

	SamplingParameters():
		nr_dir(0),
		nr_psi(0),
		nr_trans(0),
		nr_oversampled_rot(0),
		nr_oversampled_trans(0),
		nr_particles(0),
		current_oversampling(0),
		current_image_size(0),
		iclass_min(0), iclass_max(0),
		idir_min(0), idir_max(0),
		ipsi_min(0), ipsi_max(0),
		itrans_min(0), itrans_max(0),
		current_img()
	{};
};

class Indices
{
public:
	int fineIdx,
	coarseIdx,
	iclass,
	idir,
	ipsi,
	itrans,
	ioverrot,
	iovertrans;

	Indices():
		fineIdx(0),
		coarseIdx(0),
		iclass(0),
		idir(0),
		ipsi(0),
		itrans(0),
		ioverrot(0),
		iovertrans(0)
	{};

	void fineIndexToFineIndices(SamplingParameters sp) // converts an "ihidden_over" (finely sampled) index to partial indices (and coarse index)
	{
		int oversamples = sp.nr_oversampled_rot*sp.nr_oversampled_trans;
		int t_idx = fineIdx;
		iclass = floor( t_idx / ( sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples ));
		t_idx   -= iclass     * ( sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples );
		idir   = floor( t_idx / ( sp.nr_psi * sp.nr_trans * oversamples ));
		t_idx   -= idir       * ( sp.nr_psi * sp.nr_trans * oversamples );
		ipsi   = floor( t_idx / ( sp.nr_trans * oversamples ));
		t_idx   -= ipsi       * ( sp.nr_trans * oversamples );
		itrans = floor( t_idx /  oversamples );
		t_idx   -= itrans     *  oversamples ;
		ioverrot = floor( t_idx / sp.nr_oversampled_trans );
		t_idx   -= ioverrot  *   sp.nr_oversampled_trans ;
		iovertrans = t_idx ;

		coarseIdx = sp.nr_trans * sp.nr_psi * idir   +   sp.nr_trans * ipsi   +   itrans;
	}

	void fineIndicesToFineIndex(SamplingParameters sp) // converts partial indices to an "ihidden_over" (finely sampled) index // FIXME Untested
	{
		int oversamples = sp.nr_oversampled_rot*sp.nr_oversampled_trans;
		int idx = 0;
		idx += iclass   * sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples;
		idx += idir     * sp.nr_psi * sp.nr_trans * oversamples;
		idx += ipsi     * sp.nr_trans * oversamples;
		idx += itrans   * oversamples;
		idx += ioverrot * sp.nr_oversampled_trans;
		idx += iovertrans;
		fineIdx = idx;
	}

	void coarseIndexToCoarseIndices(SamplingParameters sp) // converts an "ihidden" (coarsely sampled) index to coarse partial indices // FIXME Untested
	{
		int t_idx = coarseIdx;
		iclass = floor( t_idx / ( sp.nr_dir * sp.nr_psi * sp.nr_trans));
		t_idx   -= iclass     * ( sp.nr_dir * sp.nr_psi * sp.nr_trans);
		idir   = floor( t_idx / ( sp.nr_psi * sp.nr_trans ));
		t_idx   -= idir       * ( sp.nr_psi * sp.nr_trans  );
		ipsi   = floor( t_idx / ( sp.nr_trans ));
		t_idx   -= ipsi       * ( sp.nr_trans  );
		itrans = t_idx ;
		ioverrot   = 0;
		iovertrans = 0;
	}

	void coarseIndicesToCoarseIndex(SamplingParameters sp) // converts coarse partial indices to an "ihidden" (coarsely sampled) index // FIXME Untested
	{
		int idx = 0;
		idx += idir     * sp.nr_psi * sp.nr_trans;
		idx += ipsi     * sp.nr_trans;
		idx += itrans;
		coarseIdx = idx;
	}
};


class OptimisationParamters
{
public:
	unsigned metadata_offset;

	unsigned long my_ori_particle;

	std::vector<MultidimArray<Complex > > Fimgs, Fimgs_nomask, local_Fimgs_shifted, local_Fimgs_shifted_nomask;
	std::vector<MultidimArray<RFLOAT> > Fctfs, local_Fctfs, local_Minvsigma2s;
	std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
	std::vector<RFLOAT> directions_prior, psi_prior, local_sqrtXi2;
	std::vector<RFLOAT> highres_Xi2_imgs, min_diff2, avg_diff2;
	MultidimArray<bool> Mcoarse_significant;
	// And from storeWeightedSums
	std::vector<RFLOAT> sum_weight, significant_weight, max_weight;
	std::vector<Matrix1D<RFLOAT> > old_offset, prior;
	std::vector<MultidimArray<RFLOAT> > power_imgs;
	MultidimArray<XFLOAT> Mweight;
	std::vector<Indices> max_index;

	OptimisationParamters (unsigned nr_particles, unsigned long my_ori_particle):
		metadata_offset(0),
		my_ori_particle(my_ori_particle)
	{
		power_imgs.resize(nr_particles);
		highres_Xi2_imgs.resize(nr_particles);
		Fimgs.resize(nr_particles);
		Fimgs_nomask.resize(nr_particles);
		Fctfs.resize(nr_particles);
		old_offset.resize(nr_particles);
		prior.resize(nr_particles);
		max_index.resize(nr_particles);
	};
};

class ProjectionParams
{

public:
	std::vector< size_t > orientation_num; 					// the number of significant orientation for each class
	size_t orientationNumAllClasses;							// sum of the above
	std::vector< RFLOAT > rots, tilts, psis;
	std::vector< size_t > iorientclasses, iover_rots;

	// These are arrays which detial the number of entries in each class, and where each class starts.
	// NOTE: There is no information about which class each class_idx refers to, there is only
	// a distinction between different classes.
	std::vector< size_t > class_entries, class_idx;
	inline
	ProjectionParams():

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
	ProjectionParams(size_t classes):

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
	ProjectionParams(ProjectionParams &parent, size_t start, size_t end):
		rots(				parent.rots.begin() 			+start,  	parent.rots.begin() 			+end),
		tilts(				parent.tilts.begin() 			+start, 	parent.tilts.begin() 			+end),
		psis(				parent.psis.begin() 			+start,  	parent.psis.begin() 			+end),
		iorientclasses( 	parent.iorientclasses.begin() 	+start,  	parent.iorientclasses.begin() 	+end),
		iover_rots(			parent.iover_rots.begin() 		+start,  	parent.iover_rots.begin() 		+end),
		orientation_num(1),
		orientationNumAllClasses(0),
		class_entries(1,end-start),
		class_idx(1,0) // NOTE: this is NOT the class, but rather where in these partial PrjParams to start, which is @ 0.
	{};

public:
	// Appends new values into the projection parameters for later use.
	// class_idx is used as such:
	// the n:th class (beginning with 0:th)
	// begins @ element class_idx[n]
	// ends   @ element class_idx[n]+class_entries[n]

	void pushBackAll(size_t iclass, RFLOAT NEWrot,RFLOAT NEWtilt ,RFLOAT NEWpsi, size_t NEWiorientclasses,size_t NEWiover_rots)
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

template <int AccT>
class IndexedDataArrayMaskBase
{
public:
	// indexes of job partition
	//   every element in jobOrigin    is a reference to point to a position in a IndexedDataArray.weights array where that job starts RELATIVE to firstPos
	//   every element in jobExtent    specifies the number of weights for that job
	AccPtr<size_t, AccT> jobOrigin, jobExtent;

	size_t firstPos, lastPos; // positions in indexedDataArray data and index arrays to slice out
	size_t weightNum, jobNum; // number of weights and jobs this class
	

public:

	IndexedDataArrayMaskBase() :
		jobOrigin(), 	jobExtent(), 	firstPos(), 	lastPos(),
		weightNum(), 	jobNum()
	{}
		
	IndexedDataArrayMaskBase(CudaCustomAllocator *allocator) :
		jobOrigin(allocator), 	jobExtent(allocator), 	firstPos(), 	lastPos(),
		weightNum(), 	jobNum()
	{}
	
	void setNumberOfJobs(size_t newSize)
	{
		jobNum=newSize;
		jobOrigin.resizeHost(newSize);
		jobExtent.resizeHost(newSize);
	}

	void setNumberOfWeights(size_t newSize)
	{
		weightNum=newSize;
	}

	inline
	 ~IndexedDataArrayMaskBase()
	{
//		jobOrigin.free_host();
//		jobExtent.free_host();
	};
};

template <int AccT>
class IndexedDataArrayMask : public IndexedDataArrayMaskBase<AccT>
{};

template <>
class IndexedDataArrayMask<ACC_CUDA> : public IndexedDataArrayMaskBase<ACC_CUDA>
{
public:
	inline
	IndexedDataArrayMask(CudaCustomAllocator *allocator):
	IndexedDataArrayMaskBase<ACC_CUDA>(allocator)
	{};
};

template <>
class IndexedDataArrayMask<ACC_CPU> : public IndexedDataArrayMaskBase<ACC_CPU>
{
public:
	inline
	IndexedDataArrayMask(CudaCustomAllocator *allocator):
	IndexedDataArrayMaskBase<ACC_CPU>()
	{};
};


template <int AccT>
class IndexedDataArrayBase
{
public:
	//actual data
	AccPtr<XFLOAT, AccT> weights;

	// indexes with same length as data
	// -- basic indices ---------------------------------
	//     rot_id  = id of rot     = which of all POSSIBLE orientations                               this weight signifies
	//     rot_idx = index of rot  = which in the sequence of the determined significant orientations this weight signifies
	//   trans_id  = id of trans   = which of all POSSIBLE translations                               this weight signifies
	// -- special indices ---------------------------------
	//   ihidden_overs  =  mapping to MWeight-based indexing for compatibility
	AccPtr<size_t, AccT> rot_id, rot_idx, trans_idx, ihidden_overs;

public:	
	inline
	 IndexedDataArrayBase():
		weights(),
		rot_id(),
		rot_idx(),
		trans_idx(),
		ihidden_overs()
	{};
	
	inline
	 IndexedDataArrayBase(CudaCustomAllocator *allocator):
		weights(allocator),
		rot_id(allocator),
		rot_idx(allocator),
		trans_idx(allocator),
		ihidden_overs(allocator)
	{};

	// constructor which takes a parent IndexedDataArray and a mask to create a child
	inline
	 IndexedDataArrayBase(IndexedDataArrayBase<AccT> &parent, IndexedDataArrayMask<AccT> &mask):
		weights(		&(parent.weights.hPtr[mask.firstPos])		,&(parent.weights.dPtr[mask.firstPos])			,mask.weightNum),
		rot_id(			&(parent.rot_id.hPtr[mask.firstPos])		,&(parent.rot_id.dPtr[mask.firstPos])			,mask.weightNum),
		rot_idx(		&(parent.rot_idx.hPtr[mask.firstPos])		,&(parent.rot_idx.dPtr[mask.firstPos])			,mask.weightNum),
		trans_idx(		&(parent.trans_idx.hPtr[mask.firstPos])	,&(parent.trans_idx.dPtr[mask.firstPos])		,mask.weightNum),
		ihidden_overs(	&(parent.ihidden_overs.hPtr[mask.firstPos]),&(parent.ihidden_overs.dPtr[mask.firstPos])	,mask.weightNum)
	{
		weights.doFreeDevice=false;
		rot_id.doFreeDevice=false;
		rot_idx.doFreeDevice=false;
		trans_idx.doFreeDevice=false;
		ihidden_overs.doFreeDevice=false;

		weights.doFreeHost=false;
		rot_id.doFreeHost=false;
		rot_idx.doFreeHost=false;
		trans_idx.doFreeHost=false;
		ihidden_overs.doFreeHost=false;
	};
	
	inline
	 IndexedDataArrayBase(IndexedDataArrayBase<AccT> &parent, IndexedDataArrayMask<AccT> &mask, CudaCustomAllocator *allocator):
		weights(		&(parent.weights.hPtr[mask.firstPos])		,&(parent.weights.dPtr[mask.firstPos])			,mask.weightNum, allocator),
		rot_id(			&(parent.rot_id.hPtr[mask.firstPos])		,&(parent.rot_id.dPtr[mask.firstPos])			,mask.weightNum, allocator),
		rot_idx(		&(parent.rot_idx.hPtr[mask.firstPos])		,&(parent.rot_idx.dPtr[mask.firstPos])			,mask.weightNum, allocator),
		trans_idx(		&(parent.trans_idx.hPtr[mask.firstPos])	,&(parent.trans_idx.dPtr[mask.firstPos])		,mask.weightNum, allocator),
		ihidden_overs(	&(parent.ihidden_overs.hPtr[mask.firstPos]),&(parent.ihidden_overs.dPtr[mask.firstPos])	,mask.weightNum, allocator)
	{
		weights.doFreeDevice=false;
		rot_id.doFreeDevice=false;
		rot_idx.doFreeDevice=false;
		trans_idx.doFreeDevice=false;
		ihidden_overs.doFreeDevice=false;

		weights.doFreeHost=false;
		rot_id.doFreeHost=false;
		rot_idx.doFreeHost=false;
		trans_idx.doFreeHost=false;
		ihidden_overs.doFreeHost=false;
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
		weights.resizeHost(newSize);
		rot_id.resizeHost(newSize);
		rot_idx.resizeHost(newSize);
		trans_idx.resizeHost(newSize);
		ihidden_overs.resizeHost(newSize);
	}

	void host_alloc_all()
	{
		weights.hostAlloc();
		rot_id.hostAlloc();
		rot_idx.hostAlloc();
		trans_idx.hostAlloc();
		ihidden_overs.hostAlloc();
	}

	void device_alloc_all()
	{
		weights.deviceAlloc();
		rot_id.deviceAlloc();
		rot_idx.deviceAlloc();
		trans_idx.deviceAlloc();
		ihidden_overs.deviceAlloc();
	}

	void dual_alloc_all()
	{
		host_alloc_all();
		device_alloc_all();
	}
};

template <int AccT>
class IndexedDataArray: public IndexedDataArrayBase<AccT>
{};

template <>
class IndexedDataArray<ACC_CUDA> : public IndexedDataArrayBase<ACC_CUDA>
{
public:
	inline
	IndexedDataArray(CudaCustomAllocator *allocator):
	IndexedDataArrayBase<ACC_CUDA>(allocator)
	{};
	
	inline
	IndexedDataArray(IndexedDataArrayBase<ACC_CUDA> &parent, 
			IndexedDataArrayMask<ACC_CUDA> &mask, CudaCustomAllocator *allocator):
	IndexedDataArrayBase<ACC_CUDA>(parent, mask, allocator)
	{};
};

template <>
class IndexedDataArray<ACC_CPU> : public IndexedDataArrayBase<ACC_CPU>
{
public:
	inline
	IndexedDataArray(CudaCustomAllocator *allocator):
	IndexedDataArrayBase<ACC_CPU>()
	{};
	
	inline
	IndexedDataArray(IndexedDataArrayBase<ACC_CPU> &parent, 
			IndexedDataArrayMask<ACC_CPU> &mask, CudaCustomAllocator *allocator):
	IndexedDataArrayBase<ACC_CPU>(parent, mask)
	{};
};


#if 0

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
#endif // if 0
#endif
