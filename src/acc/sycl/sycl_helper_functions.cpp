#ifdef _SYCL_ENABLED

// Make sure we build for SYCL
#include "src/acc/sycl/device_stubs.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/sycl/sycl_benchmark_utils.h"
#include "src/acc/sycl/sycl_helper_functions.h"
#include "src/acc/sycl/sycl_kernels/helper.h"
#include "src/acc/sycl/sycl_kernels/diff2_impl.h"
#include "src/acc/sycl/sycl_kernels/wavg_impl.h"
#include "src/acc/sycl/sycl_kernels/BP_impl.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"

#include "src/acc/acc_helper_functions.h"
#include "src/acc/sycl/sycl_settings.h"

#include "src/acc/acc_helper_functions_impl.h"

#endif // _SYCL_ENABLED
