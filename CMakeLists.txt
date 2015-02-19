cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

project(Relion)

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

add_executable(relion_refine src/apps/refine.cpp ${rel_SRC})
