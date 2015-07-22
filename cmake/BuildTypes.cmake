# Extra flags defined on each build type (this file is all optional to include)
#
set(CUDARCH "-arch=sm_35")



# -------------------------- 
#        Debug BUILD 
# -------------------------- 
# Additional useful nvcc-flags for debugging
#
#      -keep	            Keep all intermediate files that are generated during internal compilation steps.
#     --resource-usage      how resource usage such as registers and memeory of the GPU code. This option implies 
#                           --nvlink-options=--verbose when --relocatable-device-code=true is set. Otherwise, 
#                           it implies --ptxas-options=--verbose.

# -- Compiler flags -------------------------------------------------
set(RELION_FLAGS_DEBUG "-O0" CACHE STRING "")
set(RELION_NVCC_FLAGS_DEBUG "${CUDARCH} -lineinfo -G" CACHE STRING "")
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_DEBUG  " ")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_DEBUG        "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${RELION_LINKER_FLAGS_DEBUG}" CACHE STRING "")
set(CUDA_NVCC_FLAGS_DEBUG        "${RELION_NVCC_FLAGS_DEBUG}" CACHE STRING "")

# -- Add preprocessor defintions ------------------------------------
set(RELION_DEFINITIONS_DEBUG "-DDEBUG_CUDA -DCUDA_BENCHMARK -DUSE_THRUST")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_DEFINITIONS_DEBUG}")

#message(STATUS "Set the extra flags for Debug build type")
#message(STATUS "RELION_NVCC_FLAGS_DEBUG : ${RELION_NVCC_FLAGS_DEBUG}")
#message(STATUS "CUDA_NVCC_FLAGS_DEBUG : ${CUDA_NVCC_FLAGS_DEBUG}")
#message(STATUS "CMAKE_CXX_FLAGS_DEBUG : ${CMAKE_CXX_FLAGS_DEBUG}")
#--------------------------------------------------------------------







# -------------------------- 
#       Release BUILD 
# -------------------------- 
# Additional useful nvcc-flags for optimization
#
#     --use_fast_math
#     --prec-div         This option controls single-precision floating-point division and reciprocals. 
#                        --prec-div=true enables the IEEE round-to-nearest mode and --prec-div=false enables   
#                        the fast approximation mode. --use_fast_math implies --prec-div=false.
#     --prec-sqrt        -||- sqrt
#     --fmad             This option enables (disables) the contraction of floating-point multiplies and 
#                        adds/subtracts into floating-point multiply-add operations (FMAD, FFMA, or DFMA). 
#                        --use_fast_math implies --fmad=true.
#     --restrict         Programmer assertion that all kernel pointer parameters are restrict pointers.

# -- Compiler flags -------------------------------------------------
set(RELION_FLAGS_RELEASE "" CACHE STRING "")
set(RELION_NVCC_FLAGS_RELEASE "${CUDARCH} --disable-warnings --default-stream per-thread" CACHE STRING "")
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_RELEASE  "")

# -- Append compiler and linker flags -------------------------------
#message(STATUS "CCF_RELEASE :       ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_RELEASE}"        CACHE STRING "")
set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE} ${RELION_FLAGS_RELEASE}"          CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_RELEASE}" CACHE STRING "")
set(CUDA_NVCC_FLAGS_RELEASE        "${RELION_NVCC_FLAGS_RELEASE}"                                     CACHE STRING "")

# -- Add preprocessor defintions ------------------------------------
set(RELION_DEFINITIONS_RELEASE "-DUSE_THRUST")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_DEFINITIONS_RELEASE}")

#message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
#message(STATUS "CMAKE_CXX_FLAGS_PROFILING : ${CMAKE_CXX_FLAGS_PROFILING}")
#--------------------------------------------------------------------







# ---------------------------------- 
#       Profiling BUILD 
#  (Release with profiling output) 
# ----------------------------------
# ** NOTE: this will not have overall Release perf. **

# Additional useful nvcc-flags for profiling
#
#   -pg		               gprof profiling output (needs linker flag)
#   --resource-usage       how resource usage such as registers and memeory of the GPU code. This option implies 
#                           --nvlink-options=--verbose when --relocatable-device-code=true is set. Otherwise, 
#                           it implies --ptxas-options=--verbose#

# -- Compiler flags -------------------------------------------------
set(RELION_FLAGS_PROFILING "" CACHE STRING "")
set(RELION_NVCC_FLAGS_PROFILING "${CUDARCH} --disable-warnings -lineinfo" CACHE STRING "")
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_PROFILING  "")

# -- Append compiler and linker flags -------------------------------
#message(STATUS "CCF_RELEASE :         ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_PROFILING        "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_PROFILING}"        CACHE STRING "")
set(CMAKE_C_FLAGS_PROFILING          "${CMAKE_C_FLAGS_RELEASE} ${RELION_FLAGS_PROFILING}"        CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_PROFILING}" CACHE STRING "")
set(CUDA_NVCC_FLAGS_PROFILING        "${RELION_NVCC_FLAGS_PROFILING}"                                     CACHE STRING "")

# -- Add preprocessor defintions ------------------------------------
set(RELION_DEFINITIONS_PROFILING "-DCUDA_BENCHMARK -DTIMING -DUSE_THRUST")
set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} ${RELION_DEFINITIONS_PROFILING}")

#message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
#message(STATUS "CMAKE_CXX_FLAGS_PROFILING : ${CMAKE_CXX_FLAGS_PROFILING}")
#--------------------------------------------------------------------
