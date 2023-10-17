// For the Alternate CPU version, this is essentially a copy of
// cuda_ml_optimiser.h.  What is different is that device bundles are not
// needed, both as a separate class and referenced in MlOptimiserCpu,
// which has a few different data members and methods from MlOptimiserCuda to
// support the different implementation
// Note the the CPU implementation defines the floating point precision used
// for XFLOAT using ACC_DOUBLE_PRECISION (ACC_DOUBLE_PRECISION is also used
// for the equivalent purpose throughout the code)
#ifndef CPU_ML_OPTIMISER_H_
#define CPU_ML_OPTIMISER_H_
#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/cpu/mkl_fft.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include <stack>

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ptr.h"

class MlDataBundle
{
public:
	//The CPU accelerated projector set
	std::vector< AccProjector > projectors;

	//The CPU accelerated back-projector set
	std::vector< AccBackprojector > backprojectors;

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;
	std::vector< AccProjectorPlan > coarseProjectionPlans;

	void setup(MlOptimiser *baseMLO);

	MlDataBundle() : generateProjectionPlanOnTheFly {false}
	{}

	~MlDataBundle()
	{
		projectors.clear();
		backprojectors.clear();
	}
};


class MlOptimiserCpu
{
public:
	// transformer as holder for reuse of fftw_plans
	FourierTransformer transformer;

	MklFFT transformer1;
	MklFFT transformer2;

	MlOptimiser *baseMLO;

	bool refIs3D;
	bool dataIs3D;
    bool shiftsIs3D;

	int thread_id;

	MlDataBundle *bundle;
	std::vector< int > classStreams;


#ifdef TIMING_FILES
	relion_timer timer;
#endif

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;

	MlOptimiserCpu(MlOptimiser *baseMLOptimiser, MlDataBundle *b, const char * timing_fnm) :
			baseMLO(baseMLOptimiser),
			transformer1(baseMLOptimiser->mymodel.data_dim),
			transformer2(baseMLOptimiser->mymodel.data_dim),
			refIs3D(baseMLO->mymodel.ref_dim == 3),
            dataIs3D(baseMLO->mymodel.data_dim == 3),
            shiftsIs3D(baseMLO->mymodel.data_dim == 3 || baseMLO->mydata.is_tomo),
#ifdef TIMING_FILES
			timer(timing_fnm),
#endif
			generateProjectionPlanOnTheFly(b->generateProjectionPlanOnTheFly),
			thread_id(-1),
			bundle(b),
			classStreams(0)
	{
	};
	
	void resetData();

    void expectationOneParticle(unsigned long my_ori_particle, int thread_id);
	
	void *getAllocator()	
	{
		return nullptr;
	};

	~MlOptimiserCpu()
	{}

};
#endif
