cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

if(POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif()
project(Relion VERSION "1.4.0")

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)




# ---------------------------------------------------------SET SPECIFIC BUILD TYPE--
if(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release Profiling Benchmarking.")
else()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release Profiling Benchmarking.")
endif()
message(STATUS "Set the build type to ${CMAKE_BUILD_TYPE}")




# ------------------OPTIONS WHICH ARE NEEDED TO SET BUILD-TYPES (COMPILATION FLAGS)--
# ------------------------------------------------------------------------CUDA-ARCH--
if(NOT DEFINED CUDA_ARCH)
    message(STATUS "Setting fallback CUDA_ARCH=35")
    set(CUDARCH "-arch=sm_35")
else(NOT DEFINED CUDA_ARCH)
    message(STATUS "Using provided CUDA_ARCH=${CUDA_ARCH}")
    set(CUDARCH "-arch=sm_${CUDA_ARCH}")
endif(NOT DEFINED CUDA_ARCH)

# -----------------------------------------------DOUBLE PRECISION (CUDA-CODE) OR NOT--
option(CudaDouble "CudaDouble" OFF)
if(CudaDouble)
    add_definitions(-DCUDA_DOUBLE_PRECISION)
    message(STATUS "Building cuda files with double precision.")
	set(CudaTexture FALSE)
else(CudaDouble)
    # -- Will warn when double prec is inadvertently used in CUDA-code --
    #message(STATUS "ENV CUDA_VERSION :  $ENV{CUDA_VERSION}")
    if((DEFINED ENV{CUDA_VERSION}) AND (($ENV{CUDA_VERSION} VERSION_GREATER "7.5") OR 
    ($ENV{CUDA_VERSION} VERSION_EQUAL "7.5")))
        #message(STATUS "Cuda version is >= 7.5")
        set(WARN_DBL "--ptxas-options=-warn-double-usage") # cuda>=7.5
    endif()
endif(CudaDouble)

# ----------------------------------------------------------INCLUDE ALL BUILD TYPES--
 #This *has* to be AFTER project()
include(${CMAKE_SOURCE_DIR}/cmake/BuildTypes.cmake)

# -----------------------------------------------------------------------------CUDA--
# DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html
FIND_PACKAGE(CUDA)  
if(CUDA_FOUND)
    message(STATUS "Using cuda wrapper to compile....")
else(CUDA_FOUND)
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)

# ------------------------------------FURTHER OPTIONS PERTAINING TO CUDA-COMPILATION--
# ---------------------------------------------------------------USE TEXTURES OR NOT--
option(CudaTexture "CudaTexture" ON)
if(NOT CudaTexture)
    add_definitions(-DCUDA_NO_TEXTURES)
    message(STATUS "Texture interpolation is omitted.")
endif(NOT CudaTexture)
# ------------------------------------------------------------------ALLOCATOR CHOICE--
option(CustomAlloc "CustomAlloc" ON)
if(NOT CustomAlloc)
    add_definitions(-DCUDA_NO_CUSTOM_ALLOCATION)
    message(STATUS "Custom allocator is disabled.")
endif(NOT CustomAlloc)
# -------------------------------------------------------------FORCE USE OF STL-LIBS--
option(CudaForceSTL "CudaForceSTL" OFF)
if(CudaForceSTL)
    add_definitions(-DCUDA_FORCESTL)
    message(STATUS "Building cuda files wusing stl-libs for sort, min and max.")
endif(CudaForceSTL)
# ----------------------------------------------------------------------CUFFT OR NOT--
option(CUFFT "CUFFT" OFF)
if(CUFFT)
    if (NOT CUDA_FOUND)
        message(FATAL_ERROR "CUFFT required, but cuda was not found.")
    endif(NOT CUDA_FOUND)
    add_definitions(-DUSE_CUFFT)
    message(STATUS "Using cuda-based fftw libraries and functions")
endif(CUFFT)

# ------------------------------------------------------------------------GUI OR NOT--
# Skip FLTK/X11-dependent binaries or not
option(GUI "GUI" ON)
if(GUI)
    set(BUILD_GUI_TARGETS TRUE)
else(GUI)
    message(STATUS "You want to omit GUI targets")
endif(GUI)

# -------------------------------------------------------------------------------MPI--
find_package(MPI REQUIRED)
set(MPICXX_INCLUDE_PATH "${MPI_INCLUDE_PATH}/openmpi/ompi/mpi/cxx")
include_directories("${MPI_INCLUDE_PATH}")
include_directories("${MPICXX_INCLUDE_PATH}")
message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "MPI_LIBRARIES : ${MPI_LIBRARIES}")
message(STATUS "MPICXX_INCLUDE_PATH : ${MPICXX_INCLUDE_PATH}")
#message(STATUS "MPICXX_LIBRARIES : ${MPICXX_LIBRARIES}")

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ------------------------------------------------------------------------------FFTW--
FIND_PACKAGE(FFTW REQUIRED)

# --------------------------------------------------------------------------X11/FLTK--
#--BELOW DOES NOT WORK ON K80-MACHINES YET--
FIND_PACKAGE(X11)
if(X11_FOUND)
    set(FLTK_SKIP_OPENGL TRUE) #OpenGL is not required for relion
    FIND_PACKAGE(FLTK)
    if(FLTK_FOUND AND BUILD_GUI_TARGETS)
        message(STATUS "Both X11 and FLTK found, will compile GUI targets.")
    else()
        message(STATUS "FLTK not found or intentionally skipped; omitting GUI targets.")
    endif()
else(X11_FOUND)
    set(BUILD_GUI_TARGETS FALSE)
    message(STATUS "X11 not found, omitting GUI targets.")
endif(X11_FOUND)

# -----------------------------------------------------------------RELION COMPONENTS--
ADD_SUBDIRECTORY(src/apps)

# -----------------------------------------------------------------------------TESTS--
# Include testing flag(s) as precomiler 
# definitions and include test directives
enable_testing()
include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)

# ----------------------------------------------------------PRINT OUT ALL CMAKE VARS--
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

