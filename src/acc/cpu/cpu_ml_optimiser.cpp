// A large amount of this code is direct from cuda_ml_optimizer and so could
// be shared (but possibly with difficulty since it is enough different that
// we either need a lot of #ifdefs, or a lot of macros/some other mechanism to
// abstract the differences).  The biggest differences are the type of memory
// objects used (std::vector vs. CudaGlobalPtr and CudaCustomAllocator), the
// lack of transfers to/from the device, and on-device operations (which are
// replaced by loops/function calls).
//
// CudaFFT has been replaced with lib FFTW, if RELION is configured with mix
// precision, both single and double precision FFTW are linked into RELION.
// Install fftw-static.x86_64 and fftw-static.i686 to get the libraries without
// having to pull them at build time.  Over time we hope to replace FFTW with
// MKL.
//
// All Cuda kernels in gpu_utils and gpu_utils/cuda_kernels have been converted
// to C functions
//
// Hot spot loops in the converted C functions have been vectorized with ICC
// auto-vectorization with or without #pragma. Loop layout has been modified
// to get the best performance on CPU.
//
// NOTE:  Since the GPU code was ported back to CPU there may be additional
// changes made in the CUDA code which may not have made it here.
#ifdef ALTCPU

// Make sure we build for CPU
#include "src/acc/cpu/cuda_stubs.h"

#include "src/ml_optimiser.h"

#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <vector>
#include <iostream>
#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include "src/acc/cpu/cpu_helper_functions.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/BP.h"
#include "src/acc/cpu/mkl_fft.h"
#include "src/acc/data_types.h"
#include "src/complex.h"
#include "src/helix.h"
#include <fstream>
#include "src/parallel.h"
#include <signal.h>
#include <map>

#include <tbb/parallel_for.h>
#include <tbb/queuing_mutex.h>

#include "src/acc/utilities.h"
#include "src/acc/utilities_impl.h"

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/cpu/cpu_ml_optimiser.h"
#include "src/acc/acc_helper_functions.h"

#include "src/acc/acc_ml_optimiser_impl.h"


#include <tbb/spin_mutex.h>
tbb::spin_mutex      mkl_mutex;


void MlDataBundle::setup(MlOptimiser *baseMLO)
{
	/*======================================================
				  PROJECTOR AND BACKPROJECTOR
	======================================================*/

	unsigned nr_proj = baseMLO->mymodel.PPref.size();
	unsigned nr_bproj = baseMLO->wsum_model.BPref.size();

	projectors.resize(nr_proj);
	backprojectors.resize(nr_bproj);

	//Loop over classes
	for (int imodel = 0; imodel < nr_proj; imodel++)
	{
		projectors[imodel].setMdlDim(
				baseMLO->mymodel.PPref[imodel].data.xdim,
				baseMLO->mymodel.PPref[imodel].data.ydim,
				baseMLO->mymodel.PPref[imodel].data.zdim,
				baseMLO->mymodel.PPref[imodel].data.yinit,
				baseMLO->mymodel.PPref[imodel].data.zinit,
				baseMLO->mymodel.PPref[imodel].r_max,
				baseMLO->mymodel.PPref[imodel].padding_factor);

		projectors[imodel].initMdl(baseMLO->mdlClassComplex[imodel]);

	}

	for (int imodel = 0; imodel < nr_bproj; imodel++)
	{
		backprojectors[imodel].setMdlDim(
				baseMLO->wsum_model.BPref[imodel].data.xdim,
				baseMLO->wsum_model.BPref[imodel].data.ydim,
				baseMLO->wsum_model.BPref[imodel].data.zdim,
				baseMLO->wsum_model.BPref[imodel].data.yinit,
				baseMLO->wsum_model.BPref[imodel].data.zinit,
				baseMLO->wsum_model.BPref[imodel].r_max,
				baseMLO->wsum_model.BPref[imodel].padding_factor);

		backprojectors[imodel].initMdl();
	}

	/*======================================================
						PROJECTION PLAN
	======================================================*/

	unsigned nr_classes = baseMLO->mymodel.nr_classes;
	coarseProjectionPlans.resize(nr_classes);

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (!baseMLO->do_skip_align && !baseMLO->do_skip_rotate && !baseMLO->do_auto_refine && baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			//If doing predefined projector plan at all and is this class significant
			if (baseMLO->mymodel.pdf_class[iclass] > 0.)
			{
				std::vector<int> exp_pointer_dir_nonzeroprior;
				std::vector<int> exp_pointer_psi_nonzeroprior;
				std::vector<RFLOAT> exp_directions_prior;
				std::vector<RFLOAT> exp_psi_prior;

				long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
				long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
				long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

				coarseProjectionPlans[iclass].setup(
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


void MlOptimiserCpu::resetData()
{
	transformer1.clear();
	transformer2.clear();

	classStreams.resize(baseMLO->mymodel.nr_classes, 0);
};

void MlOptimiserCpu::expectationOneParticle(unsigned long my_part_id, int thread_id)
{
	AccPtrFactory ptrFactory(AccType::accCPU);
	accDoExpectationOneParticle<MlOptimiserCpu>(this, my_part_id, thread_id, ptrFactory);
};

#endif // ALTCPU
