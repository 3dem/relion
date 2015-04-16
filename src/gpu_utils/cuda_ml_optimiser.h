#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_

#include "src/ml_optimiser.h"


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
	MultidimArray<double> Mweight;

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

class MlOptimiserCuda
{
public:
	MlOptimiser *baseMLO;

	MlOptimiserCuda(MlOptimiser *baseMLOptimiser) : baseMLO(baseMLOptimiser) {};

	void doThreadExpectationSomeParticles(unsigned thread_id);

	void getAllSquaredDifferences(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp);

	void convertAllSquaredDifferencesToWeights(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp);

	void storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp);
};

#endif
