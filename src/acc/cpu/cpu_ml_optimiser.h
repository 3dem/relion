// For the Alternate CPU version, this is essentially a copy of
// cuda_ml_optimiser.h.  What is different is that device bundles are not
// needed, both as a separate class and referenced in MlOptimiserCpu,
// which has a few different data members and methods from MlOptimiserCuda to
// support the different implementation
// Note the the CPU implementation defines the floating point precision used
// for XFLOAT using CUDA_DOUBLE_PRECISION (CUDA_DOUBLE_PRECISION is also used
// for the equivalent purpose throughout the code)
#ifndef CPU_ML_OPTIMISER_H_
#define CPU_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/cpu/cpu_projector_plan.h"
#include "src/acc/cpu/cpu_projector.h"
#include "src/acc/cpu/cpu_backprojector.h"
#include "src/acc/cpu/mkl_fft.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include <stack>

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ptr.h"

#if 0
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


class IndexedDataArrayMask
{
public:
	// indexes of job partition
	//   every element in jobOrigin    is a reference to point to a position in a IndexedDataArray.weights array where that job starts RELATIVE to firstPos
	//   every element in jobExtent    specifies the number of weights for that job
	std::vector<size_t> jobOrigin, jobExtent;

	size_t firstPos, lastPos; // positions in indexedDataArray data and index arrays to slice out
	size_t weightNum, jobNum; // number of weights and jobs this class

	inline
	 IndexedDataArrayMask():
	firstPos(),
	lastPos(),
	weightNum(),
	jobNum()
	{};

public:

	void setNumberOfJobs(size_t newSize)
	{
		jobNum=newSize;
		jobOrigin.resize(newSize);
		jobExtent.resize(newSize);
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
	XFLOAT *weights;
	size_t  weights_size;
	bool    weights_free_memory;

	// indexes with same length as data
	// -- basic indices ---------------------------------
	//     rot_id  = id of rot     = which of all POSSIBLE orientations                               this weight signifies
	//     rot_idx = index of rot  = which in the sequence of the determined significant orientations this weight signifies
	//   trans_id  = id of trans   = which of all POSSIBLE translations                               this weight signifies
	// -- special indices ---------------------------------
	//   ihidden_overs  =  mapping to MWeight-based indexing for compatibility
	std::vector<size_t> rot_id, rot_idx, trans_idx, ihidden_overs;

	inline 
	IndexedDataArray()
	{
        weights_free_memory = false;
	};

	// constructor which takes a parent IndexedDataArray and a mask to create a child
	inline
	 IndexedDataArray(IndexedDataArray &parent, IndexedDataArrayMask &mask)
	{
              weights = parent.weights + mask.firstPos;
              weights_size = mask.weightNum;
              weights_free_memory = false;
              rot_id.resize(mask.weightNum);
              rot_idx.resize(mask.weightNum);
              trans_idx.resize(mask.weightNum);
              ihidden_overs.resize(mask.weightNum);              
              memcpy(&rot_id[0], &(parent.rot_id[mask.firstPos]), mask.weightNum * sizeof(size_t));
              memcpy(&rot_idx[0], &(parent.rot_idx[mask.firstPos]), mask.weightNum * sizeof(size_t));     
              memcpy(&trans_idx[0], &(parent.trans_idx[mask.firstPos]), mask.weightNum * sizeof(size_t));
              memcpy(&ihidden_overs[0], &(parent.ihidden_overs[mask.firstPos]), mask.weightNum * sizeof(size_t));         
	};

    ~IndexedDataArray(){
        if( weights_free_memory )
            free(weights);
    }
    
public:

	void setDataSize(size_t newSize)
	{
		weights_size = newSize;
		rot_id.resize(newSize);
		rot_idx.resize(newSize);
		trans_idx.resize(newSize);
		ihidden_overs.resize(newSize);
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
			orientation_num(1),
		orientationNumAllClasses(0),
		class_entries(1,end-start),
		class_idx(1,0) // NOTE: this is NOT the class, but rather where in these partial PrjParams to start, which is @ 0.
	{
        rots.insert(rots.begin(), parent.rots.begin()+start, parent.rots.begin()+end);
        tilts.insert(tilts.begin(), parent.tilts.begin()+start, parent.tilts.begin()+end);
        psis.insert(psis.begin(), parent.psis.begin()+start, parent.psis.begin()+end);
        iorientclasses.insert(iorientclasses.begin(), parent.iorientclasses.begin()+start, parent.iorientclasses.begin()+end);
        iover_rots.insert(iover_rots.begin(), parent.iover_rots.begin()+start, parent.iover_rots.begin()+end);	
	}

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

#endif


class MlOptimiserCpu
{
public:
	// transformer as holder for reuse of fftw_plans
	FourierTransformer transformer;

	MklFFT transformer1;
	MklFFT transformer2;

	MlOptimiser *baseMLO;

	bool refIs3D;
	bool dataIs3D;

	unsigned failsafe_attempts;

#ifdef TIMING_FILES
	relion_timer timer;
#endif

	//The CPU accelerated projector set
	std::vector< CpuProjector > cpuProjectors;

	//The CPU accelerated back-projector set
	std::vector< CpuBackprojector > cpuBackprojectors;

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;
	
	std::vector< CpuProjectorPlan > coarseProjectionPlans;
	

	MlOptimiserCpu(MlOptimiser *baseMLOptimiser, const char * timing_fnm) :
			baseMLO(baseMLOptimiser),
			transformer1(baseMLOptimiser->mymodel.data_dim),
			transformer2(baseMLOptimiser->mymodel.data_dim),
			refIs3D(baseMLO->mymodel.ref_dim == 3),
            dataIs3D(baseMLO->mymodel.data_dim == 3),
#ifdef TIMING_FILES
			timer(timing_fnm),
#endif
			failsafe_attempts(0),
			generateProjectionPlanOnTheFly(false)
	{};
	
	void resetData();

//	void doThreadExpectationSomeParticles(int thread_id);

    void expectationOneParticle(unsigned long my_ori_particle);
    
	~MlOptimiserCpu()
	{
	}

	void setupFixedSizedObjects();

	void setupTunableSizedObjects();

};

/*
class ApplyFoo {
    float *const my_a;
public:
    void operator()( const blocked_range<size_t>& r ) const {
        float *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
           Foo(a[i]);
    }
    ApplyFoo( float a[] ) :
        my_a(a)
    {}
};
*/

class cpuThreadExpectationSomeParticles {
	MlOptimiser *const my_optimiser;
public:
	void operator()( const tbb::blocked_range<size_t>& r ) const {
		MlOptimiser *mloptimiser = my_optimiser;
		MlOptimiser::CpuOptimiserType::reference ref = mloptimiser->tbbCpuOptimiser.local();
		MlOptimiserCpu *cpuOptimiser = (MlOptimiserCpu *)ref;
		if(cpuOptimiser == NULL) 
		{           
			 cpuOptimiser = new MlOptimiserCpu(mloptimiser, "cpu_optimiser");
			 cpuOptimiser->resetData();
			 cpuOptimiser->setupFixedSizedObjects();
			 cpuOptimiser->setupTunableSizedObjects();
			 ref = cpuOptimiser;
        }
		for( size_t i=r.begin(); i!=r.end(); ++i ) 
		{
			cpuOptimiser->expectationOneParticle(i);
		}
	}
	cpuThreadExpectationSomeParticles( MlOptimiser *optimiser ) :
		my_optimiser(optimiser)
	{}
};
#endif
