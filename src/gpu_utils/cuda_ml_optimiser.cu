#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/complex.h"

class CudaComplex
{
public:
	double real, imag;

	inline
	__device__ __host__ CudaComplex(): real(), imag() {};
	inline
	__device__ __host__ CudaComplex(double real, double imag): real(real), imag(imag) {};
};

class CudaImages
{
public:
	long unsigned x,y,xy,num,max_num;
	CudaComplex* start;

	inline
	__device__ __host__ CudaImages(long unsigned x, long unsigned y, long unsigned max_num):
			x(x), y(y), num(0), max_num(max_num), xy(x*y), start(new CudaComplex[xy*max_num])
	{};

	inline
	__device__ __host__ CudaComplex* current() { return start + (num*xy); };

	inline
	__device__ __host__ void increment() { num++; };

	inline
	__device__ __host__ CudaComplex* operator [](long unsigned i) { return start + (i*xy); };

	inline
	__device__ __host__ long unsigned alloc_size() { return num*xy; };

	inline
	__device__ __host__ ~CudaImages() { delete[] start; }
};

__global__ void kernel_diff2(CudaImages *ref, CudaImages *img, CudaComplex *Minvsigma2, double *partial_sums)
{
	//Dummy for now
}

void MlOptimiserCUDA::getAllSquaredDifferences(
		long int my_ori_particle, int exp_current_image_size,
		int exp_ipass, int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		std::vector<double> &exp_min_diff2,
		std::vector<double> &exp_highres_Xi2_imgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<double> > &exp_Fctfs,
		MultidimArray<double> &exp_Mweight,
		MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
		std::vector<MultidimArray<Complex> > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
		std::vector<MultidimArray<double> > &exp_local_Fctfs,
		std::vector<double> &exp_local_sqrtXi2)
{
	// Initialise min_diff and exp_Mweight for this pass
	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	exp_Mweight.resize(exp_nr_particles, mymodel.nr_classes * exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans);
	exp_Mweight.initConstant(-999.);
	if (exp_ipass==0)
		exp_Mcoarse_significant.clear();

	exp_min_diff2.clear();
	exp_min_diff2.resize(exp_nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	precalculateShiftedImagesCtfsAndInvSigma2s(false, my_ori_particle, exp_current_image_size, exp_current_oversampling,
			exp_itrans_min, exp_itrans_max, exp_Fimgs, dummy, exp_Fctfs, exp_local_Fimgs_shifted, dummy,
			exp_local_Fctfs, exp_local_sqrtXi2, exp_local_Minvsigma2s);

	// Loop only from exp_iclass_min to exp_iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		if (mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_otfshift;
			double *Minvsigma2;
			Matrix2D<double> A;

			CudaImages Frefs(exp_local_Minvsigma2s[0].xdim, exp_local_Minvsigma2s[0].ydim,
					(exp_idir_max - exp_idir_min + 1) * (exp_ipsi_max - exp_ipsi_min + 1) * exp_nr_oversampled_rot);

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses(Frefs.max_num), iover_rots(Frefs.max_num);
			long unsigned orientation_num(0);

			/*=======================================================================================
			                           Generate Reference Projections
			=========================================================================================*/

			Fref.resize(exp_local_Minvsigma2s[0]);
			Frefctf.resize(exp_local_Minvsigma2s[0]);
			if (do_shifts_onthefly)
				Fimg_otfshift.resize(Frefctf);

			std::cerr << "This is calculation of the reference projections." << std::endl;

			for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
			{
				for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
				{
					long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

					// Get prior for this direction and skip calculation if prior==0
					double pdf_orientation;
					if (do_skip_align || do_skip_rotate)
					{
						pdf_orientation = mymodel.pdf_class[exp_iclass];
					}
					else if (mymodel.orientational_prior_mode == NOPRIOR)
					{
						pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
					}
					else
					{
						pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
					}
					// In the first pass, always proceed
					// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
					// if so, proceed with projecting the reference in that direction
					bool do_proceed = (exp_ipass==0) ? true :
						isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant);
					if (do_proceed && pdf_orientation > 0.)
					{
						// Now get the oversampled (rot, tilt, psi) triplets
						// This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
						sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
								exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
						// Loop over all oversampled orientations (only a single one in the first pass)
						for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
						{
							// Get the Euler matrix
							Euler_angles2matrix(oversampled_rot[iover_rot],
												oversampled_tilt[iover_rot],
												oversampled_psi[iover_rot], A);

							Fref.data = (Complex*) Frefs.current();

							// Project the reference map (into Fref)
							(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

							//TODO REMOVE ONCE YOU KNOW THIS IS ALLWAYS TRUE
							if (Frefs.x != Fref.xdim || Frefs.y != Fref.ydim)
								std::cerr << "!!!!!!! BAD Fref size x:" << Fref.xdim << ":" << Frefs.x << " y:" << Fref.ydim << ":" << Frefs.y << std::endl;

							Frefs.increment();
							iorientclasses.push_back(iorientclass);
							iover_rots.push_back(iover_rot);
							orientation_num ++;
						}
					}
				}
			}

			CudaImages *d_Frefs;

			cudaMalloc( (void**) &d_Frefs, Frefs.alloc_size());
			cudaMemcpy( d_Frefs, Frefs.start, Frefs.alloc_size(), cudaMemcpyHostToDevice);

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CudaImages Fimgs(Frefs.x, Frefs.y,
						orientation_num * ( exp_itrans_max - exp_itrans_min + 1) * exp_nr_oversampled_trans);

				long unsigned translation_num(0);

				for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++)
				{
					sampling.getTranslations(itrans, exp_current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (exp_current_oversampling == 0)
						{
							myAB = (Frefs.y == coarse_size) ? global_fftshifts_ab_coarse[itrans].data
									: global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * exp_nr_oversampled_trans +  iover_trans;
							myAB = (strict_highres_exp > 0.) ? global_fftshifts_ab2_coarse[iitrans].data
									: global_fftshifts_ab2_current[iitrans].data;
						}
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
						{
							double real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
							double imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;

							*(Fimgs.current() + n) = CudaComplex(real, imag);
						}
						Fimgs.increment();
						translation_num ++;
					}
				}

				/*====================================
				   Initiate Particle Related On GPU
				======================================*/

				Minvsigma2 = exp_local_Minvsigma2s[ipart].data;

				CudaImages *d_Fimgs;
				double *d_Minvsigma2;

				cudaMalloc( (void**) &d_Fimgs, Fimgs.alloc_size());
				cudaMalloc( (void**) &d_Minvsigma2, Fimgs.xy);
				cudaMemcpy( d_Fimgs, Fimgs.start, Fimgs.alloc_size(), cudaMemcpyHostToDevice);
				cudaMemcpy( d_Minvsigma2, Minvsigma2, Fimgs.xy, cudaMemcpyHostToDevice);

				CudaImages *d_diff2s;
				cudaMalloc( (void**) &d_diff2s, orientation_num*translation_num);
				cudaMemset( (void**) &d_diff2s, 0, orientation_num*translation_num); //Initiate diff2 values with zeros

				/*====================================
				    		Kernel Calls
				======================================*/

				for (unsigned iorient = 0; iorient < orientation_num; iorient ++)
				{
					for (unsigned itrans = 0; itrans < translation_num; itrans ++)
					{

					}
				}

				/*====================================
				    	   Retrieve Results
				======================================*/

				double* diff2s = new double[orientation_num*translation_num];
				cudaMemcpy( diff2s, d_diff2s, orientation_num*translation_num, cudaMemcpyDeviceToHost );

				/*====================================
				    	Write To Destination TODO
				======================================*/

				/*
				for (long int i = 0; i < ihidden_overs.size(); i++)
				{
					DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_overs[i]) = diff2s[i];

					// Keep track of minimum of all diff2, only for the last image in this series
					if (diff2s[i] < exp_min_diff2[ipart])
						exp_min_diff2[ipart] = diff2s[i];
				}
				*/

				cudaFree(d_Fimgs);

			} // end loop ipart

			cudaFree(d_Frefs);

		} // end if class significant
	} // end loop iclass
}

//void MlOptimiserCUDA::precalculateModelProjectionsCtfsAndInvSigma2s(bool do_also_unmasked,
//	    int exp_current_image_size, int exp_current_oversampling,
//		std::vector<MultidimArray<Complex > > &model,
//		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
//		std::vector<MultidimArray<double> > &exp_Fctfs,
//		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
//		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
//		std::vector<MultidimArray<double> > &exp_local_Fctfs,
//		std::vector<double> &exp_local_sqrtXi2,
//		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s)
//{
//
//}

void MlOptimiserCUDA::precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked,
		long int my_ori_particle, int exp_current_image_size, int exp_current_oversampling,
		int exp_itrans_min, int exp_itrans_max,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<double> > &exp_Fctfs,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
		std::vector<MultidimArray<double> > &exp_local_Fctfs,
		std::vector<double> &exp_local_sqrtXi2,
		std::vector<MultidimArray<double> > &exp_local_Minvsigma2s)
{

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	int nr_shifts = (do_shifts_onthefly || do_skip_align) ? exp_nr_particles : exp_nr_particles * sampling.NrTranslationalSamplings(exp_current_oversampling);
	// Don't re-do if nothing has changed....
	bool do_ctf_invsig = (exp_local_Fctfs.size() > 0) ? YSIZE(exp_local_Fctfs[0])  != exp_current_image_size : true; // size has changed
	bool do_masked_shifts = (do_ctf_invsig || nr_shifts != exp_local_Fimgs_shifted.size()); // size or nr_shifts has changed

	// Use pre-sized vectors instead of push_backs!!
	exp_local_Fimgs_shifted.resize(nr_shifts);
	if (do_also_unmasked)
		exp_local_Fimgs_shifted_nomask.resize(nr_shifts);
	exp_local_Minvsigma2s.resize(exp_nr_particles);
	exp_local_Fctfs.resize(exp_nr_particles);
	exp_local_sqrtXi2.resize(exp_nr_particles);

	MultidimArray<Complex > Fimg, Fimg_nomask;
	for (int ipart = 0, my_trans_image = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);

		if (do_masked_shifts)
			windowFourierTransform(exp_Fimgs[ipart], Fimg, exp_current_image_size);
		if (do_also_unmasked)
			windowFourierTransform(exp_Fimgs_nomask[ipart], Fimg_nomask, exp_current_image_size);

		if (do_ctf_invsig)
		{
			// Also precalculate the sqrt of the sum of all Xi2
			// Could exp_current_image_size ever be different from mymodel.current_size?
			// Probably therefore do it here rather than in getFourierTransforms
			if ((iter == 1 && do_firstiter_cc) || do_always_cc)
			{
				double sumxi2 = 0.;
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
				{
					sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
				}
				// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
				exp_local_sqrtXi2[ipart] = sqrt(sumxi2);
			}

			// Also store downsized Fctfs
			// In the second pass of the adaptive approach this will have no effect,
			// since then exp_current_image_size will be the same as the size of exp_Fctfs
			windowFourierTransform(exp_Fctfs[ipart], exp_local_Fctfs[ipart], exp_current_image_size);

			// Also prepare Minvsigma2
			if (mymodel.data_dim == 3)
				exp_local_Minvsigma2s[ipart].initZeros(ZSIZE(Fimg), YSIZE(Fimg), XSIZE(Fimg));
			else
				exp_local_Minvsigma2s[ipart].initZeros(YSIZE(Fimg), XSIZE(Fimg));

			int *myMresol = (YSIZE(Fimg) == coarse_size) ? Mresol_coarse.data : Mresol_fine.data;
			// With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Minvsigma2s[ipart])
			{
				int ires = *(myMresol + n);
				// Exclude origin (ires==0) from the Probability-calculation
				// This way we are invariant to additive factors
				if (ires > 0)
					DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], n) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
			}

		}

		if (do_shifts_onthefly)
		{
			// Store a single, down-sized version of exp_Fimgs[ipart] in exp_local_Fimgs_shifted
			if (do_masked_shifts)
				exp_local_Fimgs_shifted[ipart] = Fimg;
			if (do_also_unmasked)
				exp_local_Fimgs_shifted_nomask[ipart] = Fimg_nomask;
		}
		else
		{
			// Store all translated variants of Fimg
			for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++)
			{
				// First get the non-oversampled translations as defined by the sampling object
				std::vector<double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
				sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations_x,
						oversampled_translations_y, oversampled_translations_z);
				// Then loop over all its oversampled relatives
				for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++, my_trans_image++)
				{
					// Shift through phase-shifts in the Fourier transform
					// Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
					if (do_masked_shifts)
					{
						exp_local_Fimgs_shifted[my_trans_image].resize(Fimg);
						if (mymodel.data_dim ==2)
							shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[my_trans_image],
									tab_sin, tab_cos, (double)mymodel.ori_size,
									oversampled_translations_x[iover_trans],
									oversampled_translations_y[iover_trans]);
						else
							shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[my_trans_image],
									tab_sin, tab_cos, (double)mymodel.ori_size,
									oversampled_translations_x[iover_trans],
									oversampled_translations_y[iover_trans],
									oversampled_translations_z[iover_trans]);
					}
					if (do_also_unmasked)
					{
						exp_local_Fimgs_shifted_nomask[my_trans_image].resize(Fimg_nomask);
						if (mymodel.data_dim ==2)
							shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[my_trans_image],
								tab_sin, tab_cos, (double)mymodel.ori_size,
								oversampled_translations_x[iover_trans],
								oversampled_translations_y[iover_trans]);
						else
							shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[my_trans_image],
								tab_sin, tab_cos, (double)mymodel.ori_size,
								oversampled_translations_x[iover_trans],
								oversampled_translations_y[iover_trans],
								oversampled_translations_z[iover_trans]);
					}
				}
			}
		}
	}
}


bool MlOptimiserCUDA::isSignificantAnyParticleAnyTranslation(long int iorient, int exp_itrans_min, int exp_itrans_max, MultidimArray<bool> &exp_Mcoarse_significant)
{

	long int exp_nr_trans = exp_itrans_max - exp_itrans_min + 1;
	for (long int ipart = 0; ipart < YSIZE(exp_Mcoarse_significant); ipart++)
	{
		long int ihidden = iorient * exp_nr_trans;
		for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
		{
			if (DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden))
				return true;
		}
	}
	return false;

}
