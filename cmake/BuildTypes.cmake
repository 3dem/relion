# Extra flags defined on each build type (this file is all optional to include)

# ** Debug BUILD ***************************************************************************************

# -- Compiler flags -------------------------------------------------
set(RELION_FLAGS_DEBUG "-O0 -g" CACHE STRING "")
set(RELION_NVCC_FLAGS_DEBUG "-lineinfo -G" CACHE STRING "")
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_DEBUG  " ")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${RELION_LINKER_FLAGS_DEBUG}" CACHE STRING "")
set(CUDA_NVCC_FLAGS_DEBUG "${RELION_NVCC_FLAGS_DEBUG}" CACHE STRING "")

# -- Add preprocessor defintions ------------------------------------
# This should work, but doesn't for some reason
#set(RELION_DEFINITIONS_DEBUG "CUDA_DOUBLE_PRECISION; DEBUG_CUDA; DRELION_TESTING" CACHE STRING "")
#set(COMPILE_DEFINITIONS_DEBUG "${RELION_DEFINITIONS_DEBUG} " CACHE STRING "")

set(RELION_DEFINITIONS_DEBUG "-DCUDA_DOUBLE_PRECISION -DDEBUG_CUDA -DDRELION_TESTING" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_DEFINITIONS_DEBUG}")

#message(STATUS "Set the extra flags for Debug build type")
message(STATUS "RELION_FLAGS_DEBUG : ${RELION_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG : ${CMAKE_CXX_FLAGS_DEBUG}")
#--------------------------------------------------------------------




# ** Profiling BUILD (Release with profiling output) ***************************************************
# ** NOTE: this will not have overall Release perf. ****************************************************

# -- Compiler flags -------------------------------------------------
#
#   -pg		gprof profiling output (needs linker flag)
#
set(RELION_FLAGS_PROFILING "-pg" CACHE STRING "")
set(RELION_NVCC_FLAGS_PROFILING "-lineinfo" CACHE STRING "")
# -- Linker flags ---------------------------------------------------
set(RELION_LINKER_FLAGS_PROFILING  "-pg")

# -- Append compiler and linker flags -------------------------------
set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_RELEASE} ${RELION_FLAGS_PROFILING}" CACHE STRING "")
set(CMAKE_C_FLAGS_PROFILING "${CMAKE_C_FLAGS_RELEASE} ${RELION_FLAGS_PROFILING}" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${RELION_LINKER_FLAGS_PROFILING}" CACHE STRING "")
set(CUDA_NVCC_FLAGS_PROFILING "${RELION_NVCC_FLAGS_PROFILING}" CACHE STRING "")

# -- Add preprocessor defintions ------------------------------------
set(RELION_DEFINITIONS_PROFILING "-DCUDA_DOUBLE_PRECISION " CACHE STRING "")
set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_CXX_FLAGS_PROFILING} ${RELION_DEFINITIONS_PROFILING}")

message(STATUS "RELION_FLAGS_PROFILING : ${RELION_FLAGS_PROFILING}")
message(STATUS "CMAKE_CXX_FLAGS_PROFILING : ${CMAKE_CXX_FLAGS_PROFILING}")
#--------------------------------------------------------------------