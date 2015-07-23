cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

# Set the build type
if(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
else()
   SET(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
endif()

message(STATUS "Set the build type to ${CMAKE_BUILD_TYPE}")

option(CudaDouble "CudaDouble" OFF)
if(CudaDouble)
    add_definitions(-DCUDA_DOUBLE_PRECISION)
    message(STATUS "Building cuda files with double precision.")
endif(CudaDouble)

option(CUFFT "CUFFT" OFF)
if(CUFFT)
    add_definitions(-DUSE_CUFFT)
    message(STATUS "Using cuda-based fftw libraries and functions")
endif(CUFFT)

project(Relion)
include(${CMAKE_SOURCE_DIR}/cmake/BuildTypes.cmake) #This *has* to be AFTER project()

set(PROJECT_VERSION "1.4.0")

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Uncomment following line to build relion without thread-support
# ( this will only build binareis not relying on thread-mpi support )
#set(NOTHREAD_RELION TRUE)

# Uncomment the following line to build single precision RELION 
# ( double precison is default ) 
#set(SINGLE_RELION FALSE)

# Uncomment the following to specify the code to run tests by issuing 
# mkdir build ; cd build ; cmake .. ; make ; make test


# Skip FLTK/X11-dependent binaries or not
set(OMIT_GUI_TARGETS TRUE)

find_package(CUDA)  # DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html
if(CUDA_FOUND)
    message(STATUS "Using cuda wrapper to compile....")
else(CUDA_FOUND)
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

FIND_PACKAGE(FFTW REQUIRED)

#--BELOW DOES NOT WORK ON K80-MACHINES YET--
FIND_PACKAGE(X11)

if(X11_FOUND)
    set(FLTK_SKIP_OPENGL TRUE) #OpenGL is not required for relion
    FIND_PACKAGE(FLTK)
    if(FLTK_FOUND AND NOT OMIT_GUI_TARGETS)
        message(STATUS "Both X11 and FLTK found, will compile GUI targets.")
    else()
        message(STATUS "FLTK not found or intentionally skipped; omitting GUI targets.")
    endif()
else(X11_FOUND)
    set(BUILD_GUI_TARGETS FALSE)
    message(STATUS "X11 not found, omitting GUI targets.")
endif(X11_FOUND)

ADD_SUBDIRECTORY(src/apps)

# Include testing flag(s) as precomiler definitions and include test directives
enable_testing()
include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)

# PRINT OUT ALL CMAKE VARS
#-------------------------
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#-------------------------


