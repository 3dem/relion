#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_backprojector.h"

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif

class OptimisationParamters
{
public:
	unsigned metadata_offset;

	unsigned long my_ori_particle;

	std::vector<MultidimArray<Complex > > Fimgs, Fimgs_nomask, local_Fimgs_shifted, local_Fimgs_shifted_nomask;
	std::vector<MultidimArray<double> > Fctfs, local_Fctfs, local_Minvsigma2s;
	std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
	std::vector<double> directions_prior, psi_prior, local_sqrtXi2;
	std::vector<double> highres_Xi2_imgs, min_diff2;
	MultidimArray<bool> Mcoarse_significant;
	// And from storeWeightedSums
	std::vector<double> sum_weight, significant_weight, max_weight;
	std::vector<Matrix1D<double> > old_offset, prior;
	std::vector<MultidimArray<double> > power_imgs;
	MultidimArray<XFLOAT> Mweight;

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
	};
};

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
		itrans_min(0), itrans_max(0)
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
		idir   = floor( t_idx / ( sp.nr_psi * sp.nr_trans ));
		t_idx   -= idir       * ( sp.nr_psi * sp.nr_trans  );
		ipsi   = floor( t_idx / ( sp.nr_trans ));
		t_idx   -= ipsi       * ( sp.nr_trans  );
		itrans = t_idx ;
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

class BackprojectDataBundle
{
public:
	CudaGlobalPtr<XFLOAT> reals;
	CudaGlobalPtr<XFLOAT> imags;
	CudaGlobalPtr<XFLOAT> weights;
	CudaGlobalPtr<XFLOAT> eulers;

	BackprojectDataBundle(size_t img_data_size, size_t euler_data_size, cudaStream_t stream, CudaCustomAllocator *alloc):
		reals(img_data_size, stream, alloc),
		imags(img_data_size, stream, alloc),
		weights(img_data_size, stream, alloc),
		eulers(euler_data_size, stream, alloc)
	{};
};

class MlOptimiserCuda : OutOfMemoryHandler
{
public:

	//The CUDA accelerated projector set
	std::vector< CudaProjector > cudaProjectors;

	//The CUDA accelerated back-projector set
	std::vector< CudaBackprojector > cudaBackprojectors;
	std::vector< BackprojectDataBundle *> backprojectDataBundles;

	//Used for precalculations of projection setup
	std::vector< CudaProjectorPlan > cudaCoarseProjectionPlans;

	//Used for precalculations of projection setup
	CudaCustomAllocator *allocator;

	MlOptimiser *baseMLO;

	bool refIs3D;

	bool generateProjectionPlanOnTheFly;

	int device_id;

	MlOptimiserCuda(MlOptimiser *baseMLOptimiser, int dev_id);

	void doThreadExpectationSomeParticles(unsigned thread_id);

	void storeBpMdlData()
	{
		for (int iclass = 0; iclass < baseMLO->mymodel.nr_classes; iclass++)
		{
			cudaBackprojectors[iclass].getMdlData(
					baseMLO->wsum_model.BPref[iclass].data.data,
					baseMLO->wsum_model.BPref[iclass].weight.data
					);
		}
	}

	void clearBackprojectDataBundle()
	{
		//TODO mutex lock backprojectDataBundles
		//TODO switch to cuda event synchronization instead of stream
		for (int i = 0; i < cudaBackprojectors.size(); i ++)
			cudaBackprojectors[i].syncStream();

		for (int i = 0; i < backprojectDataBundles.size(); i ++)
			delete backprojectDataBundles[i];

		backprojectDataBundles.clear();
	}

	void handleOutOfMemory()
	{
#ifdef DEBUG_CUDA
		int spaceDiff = allocator->getFreeSpace();
#endif
		clearBackprojectDataBundle();
#ifdef DEBUG_CUDA
		spaceDiff = ( (int) allocator->getFreeSpace() ) - spaceDiff;
		printf("DEBUG_INFO: MlOptimiserCuda::handleOutOfMemory called and %dB was freed.\n", spaceDiff);
#endif
	}

	~MlOptimiserCuda()
	{
		clearBackprojectDataBundle();
		delete allocator;
	}

};

#endif
