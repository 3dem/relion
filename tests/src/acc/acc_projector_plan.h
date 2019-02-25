#ifndef ACC_PROJECTOR_PLAN_H_
#define ACC_PROJECTOR_PLAN_H_

#include <vector>
#include "src/acc/acc_ptr.h"
#include "src/healpix_sampling.h"
#include <iostream>
#include <fstream>

class AccProjectorPlan
{
public:
	AccPtr< long unsigned> iorientclasses;
	AccPtr<XFLOAT> eulers;
	long unsigned orientation_num;

	AccProjectorPlan():
		orientation_num(0)
    {};
	
	AccProjectorPlan(CudaCustomAllocator *allocator):
		iorientclasses(allocator),
		eulers(allocator),
		orientation_num(0)
	{};

	//Copy constructor
	AccProjectorPlan( const AccProjectorPlan& other ):
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
			int orientational_prior_mode,
			Matrix2D<RFLOAT> &L_,
			Matrix2D<RFLOAT> &R_);

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
			int orientational_prior_mode)
	{
		Matrix2D<RFLOAT> dummyRL;

		setup(
			sampling, directions_prior, psi_prior, pointer_dir_nonzeroprior,
			pointer_psi_nonzeroprior, Mcoarse_significant, pdf_class,
			pdf_direction, nr_dir, nr_psi, nr_oversampled_rot,
			idir_min, idir_max, ipsi_min, ipsi_max, itrans_min, itrans_max,
			current_oversampling, iclass, coarse, inverseMatrix,
			do_skip_align, do_skip_rotate, orientational_prior_mode,
			dummyRL, dummyRL);
	}

	void printTo(std::ostream &os); // print

	void clear();
};

#endif
