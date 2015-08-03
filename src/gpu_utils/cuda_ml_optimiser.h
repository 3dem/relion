#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_

#include "src/ml_optimiser.h"
#include "src/gpu_utils/cuda_device_ptr.h"
#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_backprojector.h"

#ifdef CUDA_DOUBLE_PRECISION
#define FLOAT double
#else
#define FLOAT float
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
	MultidimArray<FLOAT> Mweight;

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

class MlOptimiserCuda
{
public:

	//The CUDA accelerated projector set
	std::vector< CudaProjector > cudaProjectors;

	//The CUDA accelerated back-projector set
	std::vector< CudaBackprojector > cudaBackprojectors;

	//Used for precalculations of projection setup
	std::vector< CudaProjectorPlan > cudaCoarseProjectionPlans;

	std::vector<CudaDevicePtr<FLOAT> > onTheFlyProjectionSetup;

	std::vector<CudaDevicePtr<FLOAT> > wavg_eulers;

	std::vector<CudaDevicePtr<FLOAT> > wavgs_real;
	std::vector<CudaDevicePtr<FLOAT> > wavgs_imag;
	std::vector<CudaDevicePtr<FLOAT> > wavgs_weight;

	std::vector<CudaDevicePtr<FLOAT> > Fimgs_real;
	std::vector<CudaDevicePtr<FLOAT> > Fimgs_imag;
	std::vector<CudaDevicePtr<FLOAT> > Fimgs_nomask_real;
	std::vector<CudaDevicePtr<FLOAT> > Fimgs_nomask_imag;

	std::vector<CudaDevicePtr<FLOAT> > ctfs;

	std::vector<CudaDevicePtr<FLOAT> > sorted_weights;

	std::vector<CudaDevicePtr<FLOAT> > Minvsigma2s;

	std::vector<CudaDevicePtr<FLOAT> > wdiff2s_parts;

	std::vector<CudaDevicePtr<FLOAT> > bp_eulers;

	std::vector<CudaDevicePtr<FLOAT> >     oo_otrans_x;
	std::vector<CudaDevicePtr<FLOAT> >     oo_otrans_y;
	std::vector<CudaDevicePtr<FLOAT> > myp_oo_otrans_x2y2z2;

	std::vector<CudaDevicePtr<FLOAT> >                      p_weights;
	std::vector<CudaDevicePtr<FLOAT> > p_thr_wsum_prior_offsetx_class;
	std::vector<CudaDevicePtr<FLOAT> > p_thr_wsum_prior_offsety_class;
	std::vector<CudaDevicePtr<FLOAT> >       p_thr_wsum_sigma2_offset;

	MlOptimiser *baseMLO;

	bool generateProjectionPlanOnTheFly;

	MlOptimiserCuda(MlOptimiser *baseMLOptimiser) : baseMLO(baseMLOptimiser)
	{
		unsigned nr_classes = baseMLOptimiser->mymodel.nr_classes;

		/*======================================================
                         DEVICE MEM OBJ SETUP
		======================================================*/

		wavg_eulers.resize(nr_classes);

		wavgs_real.resize(nr_classes);
		wavgs_imag.resize(nr_classes);
		wavgs_weight.resize(nr_classes);

		Fimgs_real.resize(nr_classes);
		Fimgs_imag.resize(nr_classes);
		Fimgs_nomask_real.resize(nr_classes);
		Fimgs_nomask_imag.resize(nr_classes);

		ctfs.resize(nr_classes);

		sorted_weights.resize(nr_classes);

		Minvsigma2s.resize(nr_classes);

		wdiff2s_parts.resize(nr_classes);

		bp_eulers.resize(nr_classes);

		oo_otrans_x.resize(nr_classes);
		oo_otrans_y.resize(nr_classes);
		myp_oo_otrans_x2y2z2.resize(nr_classes);

		p_weights.resize(nr_classes);
		p_thr_wsum_prior_offsetx_class.resize(nr_classes);
		p_thr_wsum_prior_offsety_class.resize(nr_classes);
		p_thr_wsum_sigma2_offset.resize(nr_classes);

		/*======================================================
            PROJECTOR, PROJECTOR PLAN AND BACKPROJECTOR SETUP
		======================================================*/

		cudaProjectors.resize(nr_classes);
		cudaBackprojectors.resize(nr_classes);

		//Can we pre-generate projector plan and corresponding euler matrices for all particles
		if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR)
			generateProjectionPlanOnTheFly = true;
		else
		{
			generateProjectionPlanOnTheFly = false;
			cudaCoarseProjectionPlans.resize(nr_classes);
		}

		//Loop over classes
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			cudaProjectors[iclass].setMdlDim(
					baseMLO->mymodel.PPref[iclass].data.xdim,
					baseMLO->mymodel.PPref[iclass].data.ydim,
					baseMLO->mymodel.PPref[iclass].data.zdim,
					baseMLO->mymodel.PPref[iclass].data.yinit,
					baseMLO->mymodel.PPref[iclass].data.zinit,
					baseMLO->mymodel.PPref[iclass].r_max,
					baseMLO->mymodel.PPref[iclass].padding_factor);

			cudaProjectors[iclass].setMdlData(baseMLO->mymodel.PPref[iclass].data.data);

			cudaBackprojectors[iclass].setMdlDim(
					baseMLO->wsum_model.BPref[iclass].data.xdim,
					baseMLO->wsum_model.BPref[iclass].data.ydim,
					baseMLO->wsum_model.BPref[iclass].data.zdim,
					baseMLO->wsum_model.BPref[iclass].data.yinit,
					baseMLO->wsum_model.BPref[iclass].data.zinit,
					baseMLO->wsum_model.BPref[iclass].r_max,
					baseMLO->wsum_model.BPref[iclass].padding_factor);

			cudaBackprojectors[iclass].initMdl();

			//If doing predefined projector plan at all and is this class significant
			if (!generateProjectionPlanOnTheFly && baseMLO->mymodel.pdf_class[iclass] > 0.)
			{
				std::vector<int> exp_pointer_dir_nonzeroprior;
				std::vector<int> exp_pointer_psi_nonzeroprior;
				std::vector<double> exp_directions_prior;
				std::vector<double> exp_psi_prior;

				long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
				long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
				long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

				cudaCoarseProjectionPlans[iclass].setup(
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

};

#endif
