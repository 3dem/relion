/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#include "src/acc/acc_projector_plan.h"
#include "src/time.h"
#include <hip/hip_runtime.h>

#ifdef _HIP_ENABLED
    #ifdef HIP_FORCESTL
        #include "src/acc/hip/hip_utils_stl.h"
    #else
        #include "src/acc/hip/hip_utils_cub.h"
    #endif
#endif

#include "src/acc/utilities.h"

#include "src/acc/acc_projector_plan_impl.h"
