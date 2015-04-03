cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

project(Relion)
set(PROJECT_VERSION "1.4.0")

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Uncomment following line to build relion without thread-support
# ( this will only build binareis not relying on thread-mpi support )
#set(NOTHREAD_RELION TRUE)

# Uncomment the following line to build single precision RELION 
# ( double precison is default ) 
#set(SINGLE_RELION TRUE)

# Uncomment the following to specify the code to run tests by issuing 
# mkdir build ; cd build ; cmake .. ; make ; make test
set(RELION_TESTING TRUE)
if(RELION_TESTING)
   enable_testing()
   add_definitions(-DRELION_TESTING)
endif()

find_package(CUDA)
if(CUDA_FOUND)
    message(STATUS "Using cuda wrapper to compile....")
else(CUDA_FOUND)
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)

# DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

FIND_PACKAGE(FFTW)

#--BELOW DOES NOT WORK ON K80-MACHINES YET--
#FIND_PACKAGE(X11)
#set(FLTK_SKIP_OPENGL TRUE) #Not required
#FIND_PACKAGE(FLTK)


ADD_SUBDIRECTORY(src/apps)

# Include testing flag(s) as precomiler definitions and include test directives
if(RELION_TESTING)
  include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)
endif()



# PRINT OUT ALL CMAKE VARS
#-------------------------
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#-------------------------


