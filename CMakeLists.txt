
cmake_minimum_required(VERSION 2.8.8)
# When we require cmake >= 2.8.12, it will provide
# CMAKE_MINIMUM_REQUIRED_VERSION automatically, but in the meantime we
# need to set a variable, and it must have a different name.
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

# CMake modules/macros are in a subdirectory to keep this file cleaner
# This needs to be set before project() in order to pick up toolchain files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Platform)

project(Relion)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# PROJECT_VERSION should have the following structure:
# VERSION-dev[-SUFFIX] where the VERSION should have the for: vMajor.vMinor.vPatch
#
# The "-dev" suffix is important to keep because it makes possible to distinguish
# between a build from official release and a build from git release branch on a
# machine with no git.
#
# NOTE: when releasing the "-dev" suffix needs to be stripped off!
# REGRESSIONTEST_VERSION and REGRESSIONTEST_BRANCH should always be
# defined.
set(PROJECT_VERSION "1.4.0")

enable_language(CXX)

include_directories("src")
include_directories("src/apps")

include_directories("include")
include_directories("include/FL")
include_directories("include/relion-1.4")

include_directories("external")
include_directories("scripts")

add_executable(relion src/apps/maingui.cpp)

install(TARGETS relion DESTINATION bin)