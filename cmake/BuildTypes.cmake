# Extra flags defined on each build type (this file is all optional to include)

#--------------Debug BUILD-------------------------------------------
# -- Compiler flags -------------------------------------------------
#
#   -pg		gprof profiling output (needs linker flag)
#

set(RELION_FLAGS_DEBUG "-O0 -g -pg" CACHE STRING "")
set(RELION_LINKER_FLAGS_DEBUG  "-pg")

# -- Linker flags ---------------------------------------------------
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${RELION_FLAGS_DEBUG}" CACHE STRING "")
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} ${RELION_LINKER_FLAGS_DEBUG}" CACHE STRING "")

# -- Preprocessor defintions ----------------------------------------

# This should work, but doesn't for some reason
#set(RELION_DEFINITIONS_DEBUG "CUDA_DOUBLE_PRECISION; DEBUG_CUDA; DRELION_TESTING" CACHE STRING "")
#set(COMPILE_DEFINITIONS_DEBUG "${COMPILE_DEFINITIONS_DEBUG} " CACHE STRING "")
#set_property(GLOBAL PROPERTY COMPILE_DEFINITIONS_DEBUG "${RELION_DEFINITIONS_DEBUG}")

set(RELION_DEFINITIONS_DEBUG "-DCUDA_DOUBLE_PRECISION -DDEBUG_CUDA -DDRELION_TESTING" CACHE STRING "")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${RELION_DEFINITIONS_DEBUG}")

#message(STATUS "Set the extra flags for Debug build type")
message(STATUS "RELION_FLAGS_DEBUG : ${RELION_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG : ${CMAKE_CXX_FLAGS_DEBUG}")
#--------------------------------------------------------------------
