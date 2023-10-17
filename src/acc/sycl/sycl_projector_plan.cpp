#ifdef _SYCL_ENABLED

// Make sure we build for SYCL
#include "src/acc/sycl/device_stubs.h"

#include "src/acc/settings.h"
#include "src/time.h"
#include "src/ml_optimiser.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/sycl/sycl_helper_functions.h"
#include "src/acc/data_types.h"
#include "src/acc/utilities.h"
#include "src/acc/acc_projector_plan.h"

#include "src/acc/acc_projector_plan_impl.h"

#endif
