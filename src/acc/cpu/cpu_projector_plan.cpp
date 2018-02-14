#ifdef ALTCPU

// Make sure we build for CPU
#include "src/acc/cpu/cuda_stubs.h"

#include "src/acc/settings.h"
#include "src/time.h"
#include "src/ml_optimiser.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/cpu/cpu_helper_functions.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/data_types.h"
#include "src/acc/utilities.h"
#include "src/acc/acc_projector_plan.h"

#include "src/acc/acc_projector_plan_impl.h"

#endif