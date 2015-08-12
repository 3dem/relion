#ifndef CUDA_PROJECTOR_PLAN_H_
#define CUDA_PROJECTOR_PLAN_H_

#include <vector>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_device_ptr.h"
#include "src/healpix_sampling.h"
#include <iostream>
#include <fstream>

class CudaProjectorPlan
{
public:
	std::vector< long unsigned > iorientclasses, iover_rots;
	long unsigned orientation_num;
	CudaDevicePtr<XFLOAT> *d_eulers;
	bool free_device;

	//Copy constructor
	CudaProjectorPlan( const CudaProjectorPlan& other ):
		iorientclasses(other.iorientclasses),
		iover_rots(other.iover_rots),
		orientation_num(other.orientation_num),
		d_eulers(other.d_eulers),
		free_device(false)
	{};

	CudaProjectorPlan():
		orientation_num(0), d_eulers(0), free_device(false)
	{
		iorientclasses.reserve(0);
		iover_rots.reserve(0);
	};

	CudaProjectorPlan(CudaDevicePtr<XFLOAT> *d_eulers):
		orientation_num(0), d_eulers(d_eulers), free_device(false)
	{
		iorientclasses.reserve(0);
		iover_rots.reserve(0);
	};

	void setup(
			HealpixSampling &sampling,
			std::vector<double> &directions_prior,
			std::vector<double> &psi_prior,
			std::vector<int> &pointer_dir_nonzeroprior,
			std::vector<int> &pointer_psi_nonzeroprior,
			MultidimArray<bool> *Mcoarse_significant,
			std::vector<double > &pdf_class,
			std::vector<MultidimArray<double> > &pdf_direction,
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

	~CudaProjectorPlan();
};

#endif
