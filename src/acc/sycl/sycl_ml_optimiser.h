// For the SYCL version, this is essentially a mix of
// cuda_ml_optimiser.h and cpu_ml_optimiser.h.
// Note the the SYCL implementation defines the floating point precision used
// for XFLOAT using ACC_DOUBLE_PRECISION (ACC_DOUBLE_PRECISION is also used
// for the equivalent purpose throughout the code)
#ifndef SYCL_ML_OPTIMISER_H_
#define SYCL_ML_OPTIMISER_H_

#include <stack>
#include <vector>
#include <tuple>
#include <typeinfo>

#include "src/mpi.h"
#include "src/ml_optimiser.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/sycl/mkl_fft.h"
#include "src/acc/sycl/sycl_benchmark_utils.h"

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ptr.h"

#include "src/acc/sycl/sycl_virtual_dev.h"

class MlSyclDataBundle
{
public:
	//The SYCL accelerated projector set
	std::vector< AccProjector > projectors;

	//The SYCL accelerated back-projector set
	std::vector< AccBackprojector > backprojectors;

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;
	std::vector< AccProjectorPlan > coarseProjectionPlans;

	void setup(MlOptimiser *baseMLO);
	void syncAllBackprojects()	{ _devAcc->waitAll(); }
	virtualSYCL* getSyclDevice()	{ return _devAcc; }

	MlSyclDataBundle(virtualSYCL *dev);
	~MlSyclDataBundle();

private:
	virtualSYCL *_devAcc;
};

class MlOptimiserSYCL
{
public:
	MlOptimiser *baseMLO;
	MlSyclDataBundle *bundle;

	// transformer as holder for reuse of fftw_plans
	FourierTransformer transformer;

	MklFFT transformer1;
	MklFFT transformer2;

	bool refIs3D;
	bool dataIs3D;
	bool shiftsIs3D;

	int threadID;

	std::vector<virtualSYCL*> classStreams;

	static void checkDevices();
	static std::vector<virtualSYCL*> getDevices(const syclDeviceType select, const std::tuple<bool,bool,bool> syclOpt, const syclBackendType BE = syclBackendType::levelZero, const bool verbose = true);

	virtualSYCL* getSyclDevice()	{ return _devAcc; }
	bool useStream() const	{ return _useStream; }

	//Used for precalculations of projection setup
	bool generateProjectionPlanOnTheFly;

	void setupDevice();

	void resetData();

	void expectationOneParticle(unsigned long my_part_id, const int thread_id);
	void doThreadExpectationSomeParticles(const int thread_id);

	void* getAllocator()
	{
		return nullptr;
	};

	MlOptimiserSYCL(MlOptimiser *baseMLOptimiser, MlSyclDataBundle *b, const bool isStream, const char *timing_fnm);

	~MlOptimiserSYCL();

private:
	virtualSYCL *_devAcc;
	bool _useStream;
};
#endif
