#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/time.h"
#include <cuda_runtime.h>

//#define PP_TIMING
#ifdef PP_TIMING
    Timer timer;
	int TIMING_TOP = timer.setNew("setup");
	int TIMING_SAMPLING = timer.setNew(" sampling");
	int TIMING_PRIOR = timer.setNew("  prior");
	int TIMING_PROC_CALC = timer.setNew("  procCalc");
	int TIMING_PROC = timer.setNew("  proc");
	int TIMING_EULERS = timer.setNew(" eulers");
#define TIMING_TIC(id) timer.tic(id)
#define TIMING_TOC(id) timer.toc(id)
#else
#define TIMING_TIC(id)
#define TIMING_TOC(id)
#endif

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

	std::vector< RFLOAT > rots, tilts, psis;
	std::vector< RFLOAT > oversampled_rot, oversampled_tilt, oversampled_psi;

	rots.reserve(nr_dir * nr_psi * nr_oversampled_rot);
	tilts.reserve(nr_dir * nr_psi * nr_oversampled_rot);
	psis.reserve(nr_dir * nr_psi * nr_oversampled_rot);

	if (iorientclasses.getSize() < nr_dir * nr_psi * nr_oversampled_rot)
	{
		iorientclasses.free_if_set();
		iorientclasses.setSize(nr_dir * nr_psi * nr_oversampled_rot);
		iorientclasses.host_alloc();
		iorientclasses.device_alloc();
	}

	orientation_num = 0;


	TIMING_TIC(TIMING_SAMPLING);

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
				sampling.getOrientations(idir, ipsi, current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
						pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior);

				// Loop over all oversampled orientations (only a single one in the first pass)
				for (long int iover_rot = 0; iover_rot < nr_oversampled_rot; iover_rot++, ipart++)
				{
					rots.push_back(oversampled_rot[iover_rot]);
					tilts.push_back(oversampled_tilt[iover_rot]);
					psis.push_back(oversampled_psi[iover_rot]);
					iorientclasses[orientation_num] = iorientclass;

					orientation_num ++;
				}
			}
			TIMING_TOC(TIMING_PROC);
		}
	}
	TIMING_TOC(TIMING_SAMPLING);

	iorientclasses.cp_to_device();

	RFLOAT alpha(.0), beta(.0), gamma(.0);
	RFLOAT ca(.0), sa(.0), cb(.0), sb(.0), cg(.0), sg(.0);
	RFLOAT cc(.0), cs(.0), sc(.0), ss(.0);

	if (eulers.getSize() < orientation_num * 9)
	{
		eulers.free_if_set();
		eulers.setSize(orientation_num * 9);
		eulers.host_alloc();
		eulers.device_alloc();
	}

	TIMING_TIC(TIMING_EULERS);

	for (long int i = 0; i < orientation_num; i++)
	{
		alpha = DEG2RAD(rots[i]);
		beta  = DEG2RAD(tilts[i]);
		gamma = DEG2RAD(psis[i]);

		sincos(alpha, &sa, &ca);
		sincos(beta,  &sb, &cb);
		sincos(gamma, &sg, &cg);

		cc = cb * ca;
		cs = cb * sa;
		sc = sb * ca;
		ss = sb * sa;

		if(inverseMatrix)
		{
			eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
			eulers[9 * i + 1] = (-sg * cc - cg * sa) ;// * padding_factor; //10
			eulers[9 * i + 2] = ( sc )               ;// * padding_factor; //20
			eulers[9 * i + 3] = ( cg * cs + sg * ca) ;// * padding_factor; //01
			eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
			eulers[9 * i + 5] = ( ss )               ;// * padding_factor; //21
			eulers[9 * i + 6] = (-cg * sb )          ;// * padding_factor; //02
			eulers[9 * i + 7] = ( sg * sb )          ;// * padding_factor; //12
			eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
		}
		else
		{
			eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
			eulers[9 * i + 1] = ( cg * cs + sg * ca) ;// * padding_factor; //01
			eulers[9 * i + 2] = (-cg * sb )          ;// * padding_factor; //02
			eulers[9 * i + 3] = (-sg * cc - cg * sa) ;// * padding_factor; //10
			eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
			eulers[9 * i + 5] = ( sg * sb )          ;// * padding_factor; //12
			eulers[9 * i + 6] = ( sc )               ;// * padding_factor; //20
			eulers[9 * i + 7] = ( ss )               ;// * padding_factor; //21
			eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
		}
	}
	TIMING_TOC(TIMING_EULERS);

	eulers.cp_to_device();

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
