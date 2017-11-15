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
	std::vector< AccProjector > projectors;
	std::vector< AccBackprojector > backprojectors;
	std::vector< AccProjectorPlan > coarseProjectionPlans;

	void setup(MlOptimiser *baseMLO);

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
#ifdef TIMING_FILES
			timer(timing_fnm),
#endif
			generateProjectionPlanOnTheFly(false),
			thread_id(-1),
			bundle(b),
			classStreams(0)
	{
		//Can we pre-generate projector plan and corresponding euler matrices for all particles
		if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR)
			generateProjectionPlanOnTheFly = true;
		else
			generateProjectionPlanOnTheFly = false;
	};
	
	void resetData();

    void expectationOneParticle(unsigned long my_ori_particle, int thread_id);
	
	CudaCustomAllocator *getAllocator()	
	{
		return ((CudaCustomAllocator *)0);
	};

	~MlOptimiserCpu()
	{}

};

/*
class ApplyFoo {
    float *const my_a;
public:
    void operator()( const blocked_range<size_t>& r ) const {
        float *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i ) 
           Foo(a[i]);
    }
    ApplyFoo( float a[] ) :
        my_a(a)
    {}
};

// Called as follows:  
// tbb::parallel_for(tbb::blocked_range<size_t>(my_first_ori_particle, my_last_ori_particle+1), 
//     cpuThreadExpectationSomeParticles(this));
class cpuThreadExpectationSomeParticles {
	MlOptimiser *const my_optimiser;
public:
	void operator()( const tbb::blocked_range<size_t>& r ) const {
		MlOptimiser *mloptimiser = my_optimiser;
		MlOptimiser::CpuOptimiserType::reference ref = mloptimiser->tbbCpuOptimiser.local();
		MlOptimiserCpu *cpuOptimiser = (MlOptimiserCpu *)ref;
		if(cpuOptimiser == NULL) 
		{           
			 cpuOptimiser = new MlOptimiserCpu(mloptimiser, "cpu_optimiser");
			 cpuOptimiser->resetData();
			 cpuOptimiser->setupFixedSizedObjects();
			 cpuOptimiser->setupTunableSizedObjects();
			 ref = cpuOptimiser;
        }
		for( size_t i=r.begin(); i!=r.end(); ++i ) 
		{
			cpuOptimiser->expectationOneParticle(i);
		}
	}
	cpuThreadExpectationSomeParticles( MlOptimiser *optimiser ) :
		my_optimiser(optimiser)
	{}
};
 */
#endif
