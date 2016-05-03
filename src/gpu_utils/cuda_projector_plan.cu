#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/time.h"
#include <cuda_runtime.h>

#define PP_TIMING
#ifdef PP_TIMING
    Timer timer;
	int TIMING_TOP = timer.setNew("setup");
	int TIMING_SAMPLING = 	timer.setNew(" sampling");
	int TIMING_PRIOR = 		timer.setNew("  prior");
	int TIMING_PROC_CALC = 	timer.setNew("  procCalc");
	int TIMING_PROC = 		timer.setNew("  proc");
	int TIMING_GEN = 		timer.setNew("   genOri");
	int TIMING_PERTURB = 	timer.setNew("   perturb");
	int TIMING_EULERS = 	timer.setNew(" eulers");
#define TIMING_TIC(id) timer.tic(id)
#define TIMING_TOC(id) timer.toc(id)
#else
#define TIMING_TIC(id)
#define TIMING_TOC(id)
#endif


void getOrientations(HealpixSampling &sampling, long int idir, long int ipsi, int oversampling_order,
		std::vector<RFLOAT > &my_rot, std::vector<RFLOAT > &my_tilt, std::vector<RFLOAT > &my_psi,
		std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
		std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior)
{
	my_rot.clear();
	my_tilt.clear();
	my_psi.clear();
	long int my_idir, my_ipsi;
	if (sampling.orientational_prior_mode == NOPRIOR)
	{
		my_idir = idir;
		my_ipsi = ipsi;
	}
	else
	{
		my_idir = pointer_dir_nonzeroprior[idir];
		my_ipsi = pointer_psi_nonzeroprior[ipsi];
	}

	if (oversampling_order == 0)
	{
		my_rot.push_back(sampling.rot_angles[my_idir]);
		my_tilt.push_back(sampling.tilt_angles[my_idir]);
		my_psi.push_back(sampling.psi_angles[my_ipsi]);
	}
	else if (!sampling.is_3D)
	{
		// for 2D sampling, only push back oversampled psi rotations
		sampling.pushbackOversampledPsiAngles(my_ipsi, oversampling_order, 0., 0., my_rot, my_tilt, my_psi);
	}
	else
	{
		// Set up oversampled grid for 3D sampling
		Healpix_Base HealPixOver(oversampling_order + sampling.healpix_order, NEST);
		int fact = HealPixOver.Nside()/sampling.healpix_base.Nside();
		int x, y, face;
		RFLOAT rot, tilt;
		// Get x, y and face for the original, coarse grid
		long int ipix = sampling.directions_ipix[my_idir];
		sampling.healpix_base.nest2xyf(ipix, x, y, face);
		// Loop over the oversampled Healpix pixels on the fine grid
		for (int j = fact * y; j < fact * (y+1); ++j)
		{
			for (int i = fact * x; i < fact * (x+1); ++i)
			{
				long int overpix = HealPixOver.xyf2nest(i, j, face);
								// this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
				double zz, phi;
				HealPixOver.pix2ang_z_phi(overpix, zz, phi);
				rot = RAD2DEG(phi);
				tilt = ACOSD(zz);

				// The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
				sampling.checkDirection(rot, tilt);

				sampling.pushbackOversampledPsiAngles(my_ipsi, oversampling_order, rot, tilt, my_rot, my_tilt, my_psi);
			}
		}
	}
}

void CudaProjectorPlan::setup(
		HealpixSampling &sampling,
		std::vector<RFLOAT> &directions_prior,
		std::vector<RFLOAT> &psi_prior,
		std::vector<int> &pointer_dir_nonzeroprior,
		std::vector<int> &pointer_psi_nonzeroprior,
		MultidimArray<bool> *Mcoarse_significant,
		std::vector<RFLOAT > &pdf_class,
		std::vector<MultidimArray<RFLOAT> > &pdf_direction,
		unsigned long nr_dir,
		unsigned long nr_psi,
		unsigned long idir_min,
		unsigned long idir_max,
		unsigned long ipsi_min,
		unsigned long ipsi_max,
		unsigned long itrans_min,
		unsigned long itrans_max,
		unsigned long current_oversampling,
		unsigned long nr_oversampled_rot,
		unsigned iclass,
		bool coarse,
		bool inverseMatrix,
		bool do_skip_align,
		bool do_skip_rotate,
		int orientational_prior_mode)
{
	TIMING_TIC(TIMING_TOP);

	std::vector< RFLOAT > oversampled_rot, oversampled_tilt, oversampled_psi;

	RFLOAT alpha(.0), beta(.0), gamma(.0);
	RFLOAT ca(.0), sa(.0), cb(.0), sb(.0), cg(.0), sg(.0);
	RFLOAT cc(.0), cs(.0), sc(.0), ss(.0);

	eulers.free_if_set();
	eulers.setSize(nr_dir * nr_psi * nr_oversampled_rot * 9);
	eulers.host_alloc();

	iorientclasses.free_if_set();
	iorientclasses.setSize(nr_dir * nr_psi * nr_oversampled_rot);
	iorientclasses.host_alloc();

	orientation_num = 0;

	TIMING_TIC(TIMING_SAMPLING);

	Matrix2D<RFLOAT> R(3,3);
	RFLOAT myperturb(0.);

	if (ABS(sampling.random_perturbation) > 0.)
	{
		myperturb = sampling.random_perturbation * sampling.getAngularSampling();
		if (sampling.is_3D)
			Euler_angles2matrix(myperturb, myperturb, myperturb, R);
	}
	else
	{
		R.initConstant(1.);
	}

	for (long int idir = idir_min, iorient = 0; idir <= idir_max; idir++)
	{
		for (long int ipsi = ipsi_min, ipart = 0; ipsi <= ipsi_max; ipsi++, iorient++)
		{
			long int iorientclass = iclass * nr_dir * nr_psi + iorient;

			TIMING_TIC(TIMING_PRIOR);
			// Get prior for this direction and skip calculation if prior==0
			RFLOAT pdf_orientation;
			if (do_skip_align || do_skip_rotate)
			{
				pdf_orientation = pdf_class[iclass];
			}
			else if (orientational_prior_mode == NOPRIOR)
			{
				pdf_orientation = DIRECT_MULTIDIM_ELEM(pdf_direction[iclass], idir);
			}
			else
			{
				pdf_orientation = directions_prior[idir] * psi_prior[ipsi];
			}
			TIMING_TOC(TIMING_PRIOR);

			// In the first pass, always proceed
			// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
			// if so, proceed with projecting the reference in that direction

			bool do_proceed(false);

			TIMING_TIC(TIMING_PROC_CALC);
			if (coarse && pdf_orientation > 0.)
				do_proceed = true;
			else if (pdf_orientation > 0.)
			{
				long int nr_trans = itrans_max - itrans_min + 1;
				for (long int ipart = 0; ipart < YSIZE(*Mcoarse_significant); ipart++)
				{
					long int ihidden = iorient * nr_trans;
					for (long int itrans = itrans_min; itrans <= itrans_max; itrans++, ihidden++)
					{
						if (DIRECT_A2D_ELEM(*Mcoarse_significant, ipart, ihidden))
						{
							do_proceed = true;
							break;
						}
					}
				}
			}
			TIMING_TOC(TIMING_PROC_CALC);

			TIMING_TIC(TIMING_PROC);
			if (do_proceed)
			{
				// Now get the oversampled (rot, tilt, psi) triplets
				// This will be only the original (rot,tilt,psi) triplet in the first pass (sp.current_oversampling==0)
				TIMING_TIC(TIMING_GEN);
				getOrientations(sampling, idir, ipsi, current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior);
				TIMING_TOC(TIMING_GEN);

				// Loop over all oversampled orientations (only a single one in the first pass)
				for (long int iover_rot = 0; iover_rot < nr_oversampled_rot; iover_rot++, ipart++)
				{
					if (sampling.is_3D)
						alpha = DEG2RAD(oversampled_rot[iover_rot]);
					else
						alpha = DEG2RAD(oversampled_rot[iover_rot] + myperturb);

				    beta  = DEG2RAD(oversampled_tilt[iover_rot]);
				    gamma = DEG2RAD(oversampled_psi[iover_rot]);

					sincos(alpha, &sa, &ca);
					sincos(beta,  &sb, &cb);
					sincos(gamma, &sg, &cg);

					cc = cb * ca;
					cs = cb * sa;
					sc = sb * ca;
					ss = sb * sa;

					Matrix2D<RFLOAT> A(3,3);

					A(0,0) = ( cg * cc - sg * sa);//00
					A(1,0) = (-sg * cc - cg * sa);//10
					A(2,0) = ( sc )              ;//20
					A(0,1) = ( cg * cs + sg * ca);//01
					A(1,1) = (-sg * cs + cg * ca);//11
					A(2,1) = ( ss )              ;//21
					A(0,2) = (-cg * sb )         ;//02
					A(1,2) = ( sg * sb )         ;//12
					A(2,2) = ( cb )              ;//22

					if (sampling.is_3D)
						A = A * R;

					if(!inverseMatrix)
						A.inv();

					eulers[9 * orientation_num + 0] = A(0,0);//00
					eulers[9 * orientation_num + 1] = A(1,0);//10
					eulers[9 * orientation_num + 2] = A(2,0);//20
					eulers[9 * orientation_num + 3] = A(0,1);//01
					eulers[9 * orientation_num + 4] = A(1,1);//11
					eulers[9 * orientation_num + 5] = A(2,1);//21
					eulers[9 * orientation_num + 6] = A(0,2);//02
					eulers[9 * orientation_num + 7] = A(1,2);//12
					eulers[9 * orientation_num + 8] = A(2,2);//22

					iorientclasses[orientation_num] = iorientclass;

					orientation_num ++;
				}
			}
			TIMING_TOC(TIMING_PROC);
		}
	}
	TIMING_TOC(TIMING_SAMPLING);

	iorientclasses.setSize(orientation_num);
	iorientclasses.put_on_device();

	eulers.setSize(orientation_num * 9);
	eulers.put_on_device();

	TIMING_TOC(TIMING_TOP);
}

void CudaProjectorPlan::printTo(std::ostream &os) // print
{
	os << "orientation_num = " << orientation_num << std::endl;
	os << "iorientclasses.size = " << iorientclasses.getSize() << std::endl;
	os << std::endl << "iorientclasses\tiover_rots\teulers" << std::endl;

	for (int i = 0; i < iorientclasses.getSize(); i ++)
	{
		os << iorientclasses[i] << "\t\t" << "\t";
		for (int j = 0; j < 9; j++)
			os << eulers[i * 9 + j] << "\t";
		os << std::endl;
	}
}

void CudaProjectorPlan::clear()
{
	orientation_num = 0;
	iorientclasses.free_if_set();
	iorientclasses.setSize(0);
	eulers.free_if_set();
	eulers.setSize(0);
#ifdef PP_TIMING
	timer.printTimes(false);
#endif
}
