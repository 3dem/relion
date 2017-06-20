#ifndef CPU_PROJECTOR_PLAN_H_
#define CPU_PROJECTOR_PLAN_H_

#include <vector>
#include "src/acc/cpu/cpu_settings.h"
#include "src/healpix_sampling.h"
#include <iostream>
#include <fstream>


class CpuProjectorPlan
{
public:
	std::vector< long unsigned> iorientclasses;
	std::vector<XFLOAT> eulers;
	long unsigned orientation_num;

	CpuProjectorPlan():
		orientation_num(0)
	{};

	//Copy constructor
	CpuProjectorPlan( const CpuProjectorPlan& other ):
		iorientclasses(other.iorientclasses),
		eulers(other.eulers),
		orientation_num(other.orientation_num)
	{};

	void setup(
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
			unsigned long nr_oversampled_rot,
			unsigned long idir_min,
			unsigned long idir_max,
			unsigned long ipsi_min,
			unsigned long ipsi_max,
			unsigned long itrans_min,
			unsigned long itrans_max,
			unsigned long current_oversampling,
			unsigned iclass,
			bool coarse,
			bool inverseMatrix,
			bool do_skip_align,
			bool do_skip_rotate,
			int orientational_prior_mode);

	void printTo(std::ostream &os); // print

	void clear();
};

#endif
