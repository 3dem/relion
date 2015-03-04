#ifndef CUDA_ML_OPTIMISER_H_
#define CUDA_ML_OPTIMISER_H_

#include <pthread.h>
#include "src/ml_model.h"
#include "src/parallel.h"
#include "src/exp_model.h"
#include "src/ctf.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/healpix_sampling.h"

class MlOptimiserCUDA
{
public:


	/* Flag to indicate orientational (i.e. rotational AND translational) searches will be skipped */
	bool do_skip_align;

	/* Flag to indicate rotational searches will be skipped */
	bool do_skip_rotate;

	// Experimental metadata model
	Experiment mydata;

	// Current ML model
	MlModel mymodel;

	// HEALPix sampling object for coarse sampling
	HealpixSampling sampling;

	// Tabulated sin and cosine functions for shifts in Fourier space
	TabSine tab_sin;
	TabCosine tab_cos;

	// Calculate translated images on-the-fly
	bool do_shifts_onthefly;

	int iter;

	// Skip marginalisation in first iteration and use signal cross-product instead of Gaussian
	bool do_firstiter_cc;

	/// Always perform cross-correlation instead of marginalization
	bool do_always_cc;

	//  Use images only up to a certain resolution in the expectation step
	int coarse_size;

    // Array with pointers to the resolution of each point in a Fourier-space FFTW-like array
	MultidimArray<int> Mresol_fine, Mresol_coarse;

	// Multiplicative fdge factor for the sigma estimates
	double sigma2_fudge;

	// Flag whether to do CTF correction
	bool do_ctf_correction;

	// Flag whether current references are ctf corrected
	bool refs_are_ctf_corrected;

	// Flag whether to do group-wise intensity bfactor correction
	bool do_scale_correction;

	std::vector<MultidimArray<Complex> > global_fftshifts_ab_coarse, global_fftshifts_ab_current,
											global_fftshifts_ab2_coarse, global_fftshifts_ab2_current;

	// Strict high-res limit in the expectation step
	double strict_highres_exp;

	MlOptimiserCUDA(Experiment mydata, MlModel mymodel, HealpixSampling sampling):
		mydata(mydata), mymodel(mymodel), sampling(sampling)
	{};

	void getAllSquaredDifferences(
			long int my_ori_particle, int exp_current_image_size,
			int exp_ipass, int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			std::vector<double> &exp_min_diff2,
			std::vector<double> &exp_highres_Xi2_imgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			MultidimArray<double> &exp_Mweight,
			MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
			std::vector<MultidimArray<double> > &exp_local_Fctfs,
			std::vector<double> &exp_local_sqrtXi2
		);

	void precalculateShiftedImagesCtfsAndInvSigma2s(
			bool do_also_unmasked,
			long int my_ori_particle, int exp_current_image_size, int exp_current_oversampling,
			int exp_itrans_min, int exp_itrans_max,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
			std::vector<MultidimArray<double> > &exp_local_Fctfs,
			std::vector<double> &exp_local_sqrtXi2,
			std::vector<MultidimArray<double> > &exp_local_Minvsigma2s
		);

	bool isSignificantAnyParticleAnyTranslation(
			long int iorient,
			int exp_itrans_min,
			int exp_itrans_max,
			MultidimArray<bool> &exp_Mcoarse_significant
		);
};

#endif
