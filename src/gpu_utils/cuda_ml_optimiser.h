#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_backprojector.h"
#include "src/gpu_utils/cuda_fft.h"
#include "src/gpu_utils/cuda_benchmark_utils.h"
#include <stack>
//#include <cufft.h>

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif

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
