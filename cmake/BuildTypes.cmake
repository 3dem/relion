# Extra flags defined on each build type (this file is all optional to include)
#
# Because gcc is compliant with a float128 type, fftw has become as well. nvcc is NOT.
# So -D__INTEL_COMPILER just manages to avoid compiling float128-targets (see fftw3.h, for instance).
# Add -G to allow cuda-gdb to break inside kernels.
if(CUDA)
    set(EXTRA_NVCC_FLAGS "-D__INTEL_COMPILER --default-stream per-thread --std=c++14")
    set(RELION_NVCC_FLAGS "${CUDARCH} ${WARN_DBL} ${EXTRA_NVCC_FLAGS}" CACHE STRING "" FORCE)
elseif (HIP)
    if (${HIP_VERSION} VERSION_LESS "5.3" )
        set(EXTRA_HIPCC_FLAGS "-fgpu-default-stream=legacy -fno-gpu-rdc -munsafe-fp-atomics")
    else()
        set(EXTRA_HIPCC_FLAGS "-fno-gpu-rdc -munsafe-fp-atomics -fgpu-default-stream=per-thread")
    endif()
    set(RELION_HIPCC_FLAGS "${EXTRA_HIPCC_FLAGS}" CACHE STRING "Compiler flags for HIP" FORCE)
endif()
# message(STATUS "RELION_NVCC_FLAGS: ${RELION_NVCC_FLAGS}")
# message(STATUS "RELION_HIPCC_FLAGS: ${RELION_HIPCC_FLAGS}")
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
set(RELION_FLAGS_DEBUG "-O0 -DDEBUG" CACHE STRING "")
if(CUDA)
    set(RELION_NVCC_FLAGS_DEBUG "${RELION_NVCC_FLAGS}" CACHE STRING "")
elseif(HIP)
    set(RELION_HIPCC_FLAGS_DEBUG "${RELION_HIPCC_FLAGS} -g -ggdb" CACHE STRING "Compiler flags for HIPCC with DEBUG configuration")
endif()

# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_DEBUG  " ")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_DEBUG        "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_DEBUG          "${CMAKE_C_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${RELION_LINKER_FLAGS_DEBUG}")

if(CUDA)
    set(CUDA_NVCC_FLAGS_DEBUG    "${RELION_NVCC_FLAGS_DEBUG}" CACHE STRING "")
elseif(HIP)
    set(CMAKE_CXX_FLAGS_DEBUG    "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_HIPCC_FLAGS_DEBUG}")
endif()

# -- Add preprocessor defintions ------------------------------------
if(CUDA)
    set(RELION_DEFINITIONS_DEBUG "-DDEBUG_CUDA")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_DEFINITIONS_DEBUG}")
elseif(HIP)
    set(RELION_DEFINITIONS_DEBUG "-DDEBUG_HIP")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_DEFINITIONS_DEBUG}")
endif()

# message(STATUS "Set the extra flags for Debug build type")
# message(STATUS "RELION_NVCC_FLAGS_DEBUG : ${RELION_NVCC_FLAGS_DEBUG}")
# message(STATUS "CUDA_NVCC_FLAGS_DEBUG : ${CUDA_NVCC_FLAGS_DEBUG}")
# message(STATUS "CMAKE_CXX_FLAGS_DEBUG : ${CMAKE_CXX_FLAGS_DEBUG}")
#--------------------------------------------------------------------




# --------------------------
#        RELWITHDEBINFO BUILD
# --------------------------
set(RELION_FLAGS_RELWITHDEBINFO "-O2 -DDEBUG" CACHE STRING "")
# -- Compiler flags -------------------------------------------------
if(CUDA)
    set(RELION_NVCC_FLAGS_RELWITHDEBINFO "${RELION_NVCC_FLAGS}" CACHE STRING "")
elseif(HIP)
    set(RELION_HIPCC_FLAGS_RELWITHDEBINFO "${RELION_HIPCC_FLAGS} -g -ggdb" CACHE STRING "Compiler flags for HIPCC with RELWITHDEBINFO configuration")
endif()

# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_RELWITHDEBINFO  " ")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${RELION_FLAGS_RELWITHDEBINFO}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO          "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${RELION_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO} ${RELION_LINKER_FLAGS_RELWITHDEBINFO}")
if(CUDA)
    set(CUDA_NVCC_FLAGS_RELWITHDEBINFO        "${RELION_NVCC_FLAGS_RELWITHDEBINFO}" CACHE STRING "")
    # -- Add preprocessor defintions ------------------------------------
    set(RELION_DEFINITIONS_RELWITHDEBINFO "-DDEBUG_CUDA")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${RELION_DEFINITIONS_RELWITHDEBINFO}")
elseif(HIP)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO        "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${RELION_HIPCC_FLAGS_RELWITHDEBINFO}")
    # -- Add preprocessor defintions ------------------------------------
    set(RELION_DEFINITIONS_RELWITHDEBINFO "-DDEBUG_HIP")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${RELION_DEFINITIONS_RELWITHDEBINFO}")
endif()


# message(STATUS "Set the extra flags for RELWITHDEBINFO build type")
# message(STATUS "RELION_NVCC_FLAGS_RELWITHDEBINFO : ${RELION_NVCC_FLAGS_RELWITHDEBINFO}")
# message(STATUS "CUDA_NVCC_FLAGS_RELWITHDEBINFO : ${CUDA_NVCC_FLAGS_RELWITHDEBINFO}")
# message(STATUS "CMAKE_CXX_FLAGS_RELWITHDEBINFO : ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
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
set(RELION_FLAGS_RELEASE " " CACHE STRING "")
if(CUDA)
    set(RELION_NVCC_FLAGS_RELEASE "${RELION_NVCC_FLAGS} --disable-warnings" CACHE STRING "")
elseif(HIP)
    set(RELION_HIPCC_FLAGS_RELEASE "${RELION_HIPCC_FLAGS} -w" CACHE STRING "Compiler flags for HIPCC with RELEASE configuration")
endif()
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_RELEASE  "")

# -- Append compiler and linker flags -------------------------------
#message(STATUS "CCF_RELEASE :       ${CMAKE_CXX_FLAGS_RELEASE}")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_RELEASE}")
set(CMAKE_C_FLAGS_RELEASE          "${CMAKE_C_FLAGS_RELEASE} ${RELION_FLAGS_RELEASE}")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_RELEASE}")
if(CUDA)
    set(CUDA_NVCC_FLAGS_RELEASE    "${RELION_NVCC_FLAGS_RELEASE}" CACHE STRING "")
elseif(HIP)
    set(CMAKE_CXX_FLAGS_RELEASE    "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_HIPCC_FLAGS_RELEASE}")
endif()
# -- Add preprocessor defintions ------------------------------------
set(RELION_DEFINITIONS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_DEFINITIONS_RELEASE}")

# message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
# message(STATUS "CMAKE_CXX_FLAGS_RELEASE : ${CMAKE_CXX_FLAGS_RELEASE}")
#--------------------------------------------------------------------


if(CUDA)
    # ----------------------------------
    #       NVIDIA Profiling BUILD
    #        (Release for nvprof)
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
    set(RELION_NVCC_FLAGS_PROFILING "${RELION_NVCC_FLAGS} -lineinfo" CACHE STRING "")
    # -- Linker flags ---------------------------------------------------
    set(RELION_LINKER_FLAGS_PROFILING  "")

    # -- Append compiler and linker flags -------------------------------
    set(CMAKE_CXX_FLAGS_PROFILING        "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_PROFILING}"               CACHE STRING "")
    set(CMAKE_C_FLAGS_PROFILING          "${CMAKE_C_FLAGS_RELEASE} ${RELION_FALAGS_PROFILING}"                 CACHE STRING "")
    set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_PROFILING}" CACHE STRING "")
    set(CUDA_NVCC_FLAGS_PROFILING        "${RELION_NVCC_FLAGS_PROFILING}"                                     CACHE STRING "")

    # -- Add preprocessor defintions ------------------------------------
    set(RELION_DEFINITIONS_PROFILING "-DCUDA_PROFILING")
    set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} ${RELION_DEFINITIONS_PROFILING}")

    #message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
    #message(STATUS "CMAKE_CXX_FLAGS_PROFILING : ${CMAKE_CXX_FLAGS_PROFILING}")
    #--------------------------------------------------------------------

elseif(HIP)
    # ----------------------------------
    #       AMD Profiling BUILD
    #        (Release for rocprof)
    # ----------------------------------
    # ** NOTE: this will not have overall Release perf. **


    # -- Compiler flags -------------------------------------------------
    set(RELION_FLAGS_PROFILING "-O2" CACHE STRING "")
    set(RELION_HIPCC_FLAGS_PROFILING "${RELION_HIPCC_FLAGS} -g -ggdb" CACHE STRING "Semicolon delimited flags")
    # -- Linker flags ---------------------------------------------------
    set(RELION_LINKER_FLAGS_PROFILING  "")

    # -- Append compiler and linker flags -------------------------------
    set(CMAKE_CXX_FLAGS_PROFILING        "${RELION_FLAGS_PROFILING} ${RELION_HIPCC_FLAGS_PROFILING}")
    set(CMAKE_C_FLAGS_PROFILING          "${RELION_FLAGS_PROFILING} ${RELION_HIPCC_FLAGS_PROFILING}")
    set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_PROFILING}")

    # -- Add preprocessor defintions ------------------------------------
    set(RELION_DEFINITIONS_PROFILING "-DHIP_PROFILING -DDEBUG_HIP -DDEBUG")
    set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} ${RELION_DEFINITIONS_PROFILING}")

    # message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
    # message(STATUS "CMAKE_CXX_FLAGS_PROFILING : ${CMAKE_CXX_FLAGS_PROFILING}")
    #--------------------------------------------------------------------
endif()



# ----------------------------------
#       Benchmarking BUILD
#  (Release with profiling output)
# ----------------------------------
# -- Compiler flags -------------------------------------------------
set(RELION_FLAGS_BENCHMARKING "" CACHE STRING "")
if(CUDA)
    set(RELION_NVCC_FLAGS_BENCHMARKING "${RELION_NVCC_FLAGS} " CACHE STRING "")
elseif(HIP)
    set(RELION_HIPCC_FLAGS_BENCHMARKING "${RELION_HIPCC_FLAGS} " CACHE STRING "Semicolon delimited flags for HIPCC with Benchmarking configuration")
endif()
    # -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_BENCHMARKING  "")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_BENCHMARKING        "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_BENCHMARKING}"               CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_BENCHMARKING          "${CMAKE_C_FLAGS_RELEASE} ${RELION_FLAGS_BENCHMARKING}"                 CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_BENCHMARKING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_BENCHMARKING}" CACHE STRING "" FORCE)
if(CUDA)
    set(CUDA_NVCC_FLAGS_BENCHMARKING    "${RELION_NVCC_FLAGS_BENCHMARKING}" CACHE STRING "" FORCE)
elseif(HIP)
    set(CMAKE_CXX_FLAGS_BENCHMARKING    "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_HIPCC_FLAGS_BENCHMARKING}")
endif()
# -- Add preprocessor defintions ------------------------------------
if(CUDA)
    set(RELION_DEFINITIONS_BENCHMARKING "-DCUDA_BENCHMARK -DTIMING")
    set(CMAKE_CXX_FLAGS_BENCHMARKING "${CMAKE_CXX_FLAGS_BENCHMARKING} ${RELION_DEFINITIONS_BENCHMARKING}")
elseif(HIP)
    set(RELION_DEFINITIONS_BENCHMARKING "-DHIP_BENCHMARK -DTIMING")
    set(CMAKE_CXX_FLAGS_BENCHMARKING "${CMAKE_CXX_FLAGS_BENCHMARKING} ${RELION_DEFINITIONS_BENCHMARKING}")
endif()
#--------------------------------------------------------------------
# message(STATUS "CMAKE_CXX_FLAGS_BENCHMARKING : ${CMAKE_CXX_FLAGS_BENCHMARKING}")
