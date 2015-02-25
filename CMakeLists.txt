cmake_minimum_required(VERSION 2.8.8)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.8")

project(Relion)
set(PROJECT_VERSION "1.4.0")

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Uncomment the lines below to force local builds of fftw and/or fltk
#set(BUILD_FFTW TRUE)
#set(BUILD_FLTK TRUE)

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

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

set(REL_EXTERNAL_LIBS_TAR_DIRECTORY ${CMAKE_SOURCE_DIR}/external)
set(REL_EXTERNAL_LIBS_EXTRACT_TARGET ${CMAKE_BINARY_DIR}/external)
include(ExternalProject)


  message( STATUS " suffs: " ${CMAKE_FIND_LIBRARY_SUFFIXES})
LIST(APPEND CMAKE_FIND_LIBRARY_SUFFIXES .so.1)
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so.1 ${CMAKE_FIND_LIBRARY_SUFFIXES}")
  message( STATUS " suffs: " ${CMAKE_FIND_LIBRARY_SUFFIXES})

#Try find FFTW -------------------------------------------------------------------------------------
if(NOT BUILD_FFTW)
  find_package(FFTW) 
  message( STATUS " FFTW: " ${FFTW_LIBRARIES})
  message( STATUS " FFTW_PATH: " ${FFTW_PATH}) 
  #message( STATUS " FFTW_FOUND: " ${FFTW_FOUND}) 
endif()

#Build FFTW if not found
if(FFTW_FOUND)
   message( STATUS "Found fftw, will NOT build")
else(FFTW_FOUND)
   message( STATUS "no fftw found, WILL build")
   set(REL_FFTW3_TAR_FILE ${REL_EXTERNAL_LIBS_TAR_DIRECTORY}/fftw-3.2.2.tar.gz)
   set(REL_FFTW3_LIB_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
   set(REL_FFTW3_BUILD_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)
   externalproject_add(FFTW3
	URL ${REL_FFTW3_TAR_FILE}
	SOURCE_DIR ${REL_FFTW3_LIB_DIR}
	INSTALL_DIR ${CMAKE_BINARY_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure --enable-threads --enable-shared --prefix=<INSTALL_DIR>
	BUILD_COMMAND ${MAKE}
	LOG_CONFIGURE
	LOG_BUILD
	LOG_INSTALL
   )
endif(FFTW_FOUND)
#------------------------------------------------------------------------------------------------------

#Try to find X11------------------------------------------------------------
find_package(X11)
if(X11_FOUND)
   message( STATUS "Found X11, will NOT build")
   message( STATUS ${X11_LIBRARIES})
endif(X11_FOUND)

find_path (X11_INCLUDES Xdbe.h)
#find_library (FFTW_LIBRARIES NAMES fftw3)
#   message( STATUS " FFTW: " ${FFTW_LIBRARIES})
#find_path(FFTW_PATH fftw3.h /usr/lib/* /usr/share/* /usr/include/* )
#   message( STATUS " FFTW_PATH: " ${FFTW_PATH}) 
#----------------------------------------------------------------------------


#Try to find FLTK--------------------------------------------------------------------------------------
if(NOT BUILD_FLTK)
   find_package(FLTK)
      message( STATUS " FLTK: " ${FLTK_LIBRARIES})
   #find_package(FLTK REQUIRED NO_MODULE) 
   include(${FLTK_USE_FILE})
   find_path(FLTK_INCLUDES Fl.H)
   find_library(FLTK_LIBRARIES NAMES FL)
      message( STATUS " FLTK: " ${FLTK_LIBRARIES})
   find_path(FL_PATH Fl.H /usr/include /usr/lib/* /usr/share/* /usr/include/* )
      message( STATUS " FL_PATH: " ${FL_PATH})
endif(NOT BUILD_FLTK)

#Build FLTK if not found
if(FLTK_FOUND)
   message( STATUS "Found fltk, will NOT build")

else(FLTK_FOUND)
   message( STATUS "No fltk found, WILL build")
   set(REL_FLTK_TAR_FILE ${REL_EXTERNAL_LIBS_TAR_DIRECTORY}/fltk-1.3.0.tar.gz)
   set(REL_FLTK_LIB_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fltk)
   set(REL_FLTK_BUILD_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fltk-build)

   externalproject_add(FLTK
	URL ${REL_FLTK_TAR_FILE}
	SOURCE_DIR ${REL_FLTK_LIB_DIR}
	INSTALL_DIR ${CMAKE_BINARY_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure --enable-shared --prefix=<INSTALL_DIR>
	BUILD_COMMAND ${MAKE}
	BUILD_IN_SOURCE
	BINARY_DIR ${REL_FLTK_LIB_DIR}
	LOG_CONFIGURE
	LOG_BUILD
	LOG_INSTALL
   )
endif(FLTK_FOUND)
#------------------------------------------------------------------------------------------------------


include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c")
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h")
file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")
file(GLOB REL_LIB "${CMAKE_SOURCE_DIR}/lib/*.so" )

file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")
#	message(STATUS "ALL: ${_target}...")

add_library(relion_lib OBJECT ${REL_SRC} ${REL_SRC_H} ${REL_HP})
foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
	if(CUDA_FOUND)
	    cuda_add_executable(${_target} src/apps/${_target}.cpp $<TARGET_OBJECTS:relion_lib>) # ${REL_SRC} ${REL_SRC_H} ${REL_HP})
	else(CUDA_FOUND)
	    add_executable(${_target} src/apps/${_target}.cpp $<TARGET_OBJECTS:relion_lib>) # ${REL_SRC} ${REL_SRC_H} ${REL_HP})
	endif(CUDA_FOUND)
	add_dependencies(${_target} FFTW FLTK)
	#target_link_libraries(${_target} fltk ${REL_LIB} libXft.so libXext.so libXinerama.so ) # ${REL_SRC_H} ${REL_SRC} ${REL_LIB}) # ${BUILT_STUFF})
	target_link_libraries(${_target} ${FLTK_LIBRARIES} ${FFTW_LIBRARIES} ${X11} ${relion_lib} ) # ${REL_SRC_H} ${REL_SRC} ${REL_LIB}) # ${BUILT_STUFF})
	message(STATUS "added ${_target}...")
endforeach()

#foreach (_target ${RELION_TARGETS})
#	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
#	add_executable(${_target} src/apps/${_target}.cpp ${REL_SRC} ${REL_SRC_H} ${REL_HP})
#	add_dependencies(${_target} FFTW FLTK)
#	#target_link_libraries(${_target} fltk ${REL_LIB} libXft.so libXext.so libXinerama.so ) # ${REL_SRC_H} ${REL_SRC} ${REL_LIB}) # ${BUILT_STUFF})
#	target_link_libraries(${_target} ${FLTK_LIBRARIES} ${FFTW_LIBRARIES} ${X11} ) # ${REL_SRC_H} ${REL_SRC} ${REL_LIB}) # ${BUILT_STUFF})
#	message(STATUS "added ${_target}...")
#endforeach()

# Include testing flag(s) as precomiler definitions and include test directives
if(RELION_TESTING)
  add_definitions(-DRELION_TESTING)
  include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)
endif()



# PRINT OUT ALL CMAKE VARS
#-------------------------
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#-------------------------

