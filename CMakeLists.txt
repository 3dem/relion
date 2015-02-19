cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

project(Relion)

find_package(CUDA)
# DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html

# PRINT OUT ALL CMAKE VARS
#-------------------------
get_cmake_property(_variableNames VARIABLES)
foreach (_variableName ${_variableNames})
    message(STATUS "${_variableName}=${${_variableName}}")
endforeach()
#-------------------------

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(PROJECT_VERSION "1.4.0")

include_directories(".")
include_directories("external/fltk-1.3.0/")

file(GLOB rel_SRC
    "src/*.cpp"
)

if(CUDA_FOUND)
    cuda_add_executable(relion_refine src/apps/refine.cpp ${rel_SRC})
    message(STATUS "Using cuda wrapper to compile....")
else(CUDA_FOUND)
    # --probably don't even want to compile in that case, but for now--
    add_executable(relion_refine src/apps/refine.cpp ${rel_SRC})
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)
