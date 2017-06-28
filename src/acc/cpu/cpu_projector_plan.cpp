#ifdef ALTCPU

// Make sure we build for CPU
#undef CUDA
typedef float cudaStream_t;
typedef double CudaCustomAllocator;
typedef int dim3;
#define cudaStreamPerThread 0
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) //Do nothing
#define LAUNCH_PRIVATE_ERROR(func, status)
#define LAUNCH_HANDLE_ERROR( err )
#define DEBUG_HANDLE_ERROR( err )

#include "src/acc/settings.h"
#include "src/time.h"

#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_helper_functions.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"
#include "src/acc/acc_ptr.h"

#include "src/acc/acc_projector_plan.h"

#include "src/acc/acc_projector_plan_impl.h"

#endif