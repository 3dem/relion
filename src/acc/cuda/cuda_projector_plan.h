#ifndef CUDA_PROJECTOR_PLAN_H_
#define CUDA_PROJECTOR_PLAN_H_

#include <vector>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/healpix_sampling.h"
#include <iostream>
#include <fstream>

class CudaProjectorPlan
{
public:
	CudaGlobalPtr< long unsigned> iorientclasses;
	CudaGlobalPtr<XFLOAT> eulers;
	long unsigned orientation_num;

	CudaProjectorPlan(CudaCustomAllocator *allocator):
		iorientclasses(allocator),
		eulers(allocator),
		orientation_num(0)
	{};

	//Copy constructor
	CudaProjectorPlan( const CudaProjectorPlan& other ):
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
