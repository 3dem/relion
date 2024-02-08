cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

#if(POLICY CMP0048)
#  cmake_policy(SET CMP0048 NEW)
#endif()
project(Relion)

# Use new policy for OS X @rpath
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_definitions(-DINSTALL_LIBRARY_DIR=${CMAKE_INSTALL_PREFIX}/lib/)
add_definitions(-DSOURCE_DIR=${CMAKE_SOURCE_DIR}/src/)

# message(STATUS "INSTALL_LIBRARY_DIR set to ${CMAKE_INSTALL_PREFIX}/lib/")
# message(STATUS "SOURCE_DIR set to ${CMAKE_SOURCE_DIR}/src/")

# ------------------------------------------------------------------RPATH SETTINGS--
if(NOT APPLE)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
        SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
endif(NOT APPLE)

# ---------------------------------------------------------SET SPECIFIC BUILD TYPE--
if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "")
	string( TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER )
	
	if( ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "none" ) AND
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "release" ) AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug" ) AND
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "relwithdebinfo" ) AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "profiling" )  AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "benchmarking" ) )
	     message( FATAL_ERROR "CMAKE_BUILD_TYPE : '${CMAKE_BUILD_TYPE}' is not a valid build type. "
		"Valid options are: 'None', 'Release', 'Debug', 'RelWithDebInfo', and 'Profiling'." )
	endif()
	
	message(STATUS "BUILD TYPE set to '${CMAKE_BUILD_TYPE}'")
	SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: 'None', 'Release', 'Debug', 'RelWithDebInfo', and 'Profiling'.")

	if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug")
		set(MDT_TYPE_CHECK ON)
	endif()

else()
	SET(CMAKE_BUILD_TYPE "Release")
	message(STATUS "BUILD TYPE set to the default type:  '${CMAKE_BUILD_TYPE}'")
	string( TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER )
endif()

# ------------------OPTIONS WHICH ARE NEEDED TO SET BUILD-TYPES (COMPILATION FLAGS)--
# ------------------------------------------------------------------------CUDA-ARCH--
if(NOT DEFINED CUDA_ARCH)
    message(STATUS "Setting fallback CUDA_ARCH=50")
    set(CUDARCH "-arch=sm_50")
else(NOT DEFINED CUDA_ARCH)
    message(STATUS "Using provided CUDA_ARCH=${CUDA_ARCH}")
    set(CUDARCH "-arch=sm_${CUDA_ARCH}")
endif(NOT DEFINED CUDA_ARCH)
    
# -------------------------------------------------------------------FURTHER OPTIONS--


# CUDA on by default, so check for CPU-accelration request and possible conflicting dual-request
option(ALTCPU "Enable Accelerated CPU version" OFF)
if(ALTCPU)
	option(CUDA "Enable CUDA GPU acceleration" OFF)
	if(ALTCPU AND CUDA)
   		message(FATAL_ERROR "You cannot build with both CUDA=ON and ALTCPU=ON.  Please choose one and rerun CMAKE")
	endif()
else()
	option(CUDA "Enable CUDA GPU acceleration" ON)
endif()

option(DoublePrec_CPU "DoublePrec_CPU" ON)
option(DoublePrec_ACC "Accelerated Code use double-precision" OFF)
option(MKLFFT "Use MKL rather than FFTW for FFT" OFF)
option(CudaTexture "CudaTexture" ON)

if(ALTCPU)
    message(STATUS "ALTCPU enabled - Building CPU-accelerated version of RELION")
endif()

if(CUDA)
    message(STATUS "CUDA enabled - Building CUDA-accelerated version of RELION")
endif()

if(CUDA OR ALTCPU)
    add_definitions(-DACC_CUDA=2 -DACC_CPU=1)
endif()

# -----------------------------------------------DOUBLE PRECISION (CUDA-CODE) OR NOT--
if(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to double")    
else(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to single")
    add_definitions(-DRELION_SINGLE_PRECISION)
endif(DoublePrec_CPU)

if(DoublePrec_ACC)
    message(STATUS "Setting accelerated code precision to double")
    add_definitions(-DACC_DOUBLE_PRECISION)
	set(CudaTexture FALSE)
else(DoublePrec_ACC)
    message(STATUS "Setting accelerated code precision to single")
endif(DoublePrec_ACC)

if(MDT_TYPE_CHECK)
    message("-- MetaDataTable type check is enabled.")
    add_definitions(-DMETADATA_TABLE_TYPE_CHECK)
endif()

# ----------------------------------------------------------INCLUDE ALL BUILD TYPES--
 #This *has* to be AFTER project()
include(${CMAKE_SOURCE_DIR}/cmake/BuildTypes.cmake)

if(CMAKE_VERSION VERSION_LESS 3.10)
	FIND_PACKAGE(CUDA)
else(CMAKE_VERSION)
	# newer version of CMake use a different module name
	FIND_PACKAGE(CUDAToolkit)
endif()

if(CUDA_FOUND)
    message(STATUS "Using cuda wrapper to compile....")
    if( (NOT ${CUDA_VERSION} VERSION_LESS "7.5") AND (NOT DoublePrec_ACC) )
        message(STATUS "Cuda version is >= 7.5 and single-precision build, enable double usage warning.")
		set(WARN_DBL "--ptxas-options=-warn-double-usage") # cuda>=7.5
	elseif( ${CUDA_VERSION} VERSION_LESS "7.0")
		message(WARNING "Cuda version is less than 7.0, so relion will be compiled without GPU support.")
		set(CUDA OFF)
    endif()
    
    if(CUDA)
        add_definitions(-D_CUDA_ENABLED)
    endif()
else(CUDA_FOUND)
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)

# ------------------------------------------------------------------ALLOCATOR CHOICE--
option(CachedAlloc "CachedAlloc" ON)
if(NOT CachedAlloc)
    add_definitions(-DCUDA_NO_CUSTOM_ALLOCATION)
	message(STATUS "Cached allocation is disabled.")
endif(NOT CachedAlloc)
option(CustomAllocMemGuards "CustomAllocMemGuards" OFF)
if(CustomAllocMemGuards)
    add_definitions(-DCUSTOM_ALLOCATOR_MEMGUARD)
    message(STATUS "Abort on out of bound write.")
endif(CustomAllocMemGuards)
# -------------------------------------------------------------FORCE USE OF STL-LIBS--
option(CudaForceSTL "CudaForceSTL" OFF)
if(CudaForceSTL)
    add_definitions(-DCUDA_FORCESTL)
    message(STATUS "Building cuda files wusing stl-libs for sort, min and max.")
endif(CudaForceSTL)

# ------------------------------------------------------------------------GUI OR NOT--
# Skip FLTK/X11-dependent binaries or not
option(GUI "GUI" ON)
if(NOT GUI)
    message(STATUS "Omitting GUI targets as per your request")
endif()

# ---------------------------------------------------------------------------TBB --

option(FORCE_OWN_TBB "FORCE_OWN_TBB" OFF)

if (ALTCPU)
	
	if (FORCE_OWN_TBB)
		message(STATUS "Will ignore any potentially installed system TBB lib, as per your request.")
		include(${CMAKE_SOURCE_DIR}/cmake/BuildTBB.cmake)
		set(INSTALL_OWN_TBB 1)
	else(FORCE_OWN_TBB)
		find_package(TBB)
		if(TBB_FOUND)
			include_directories("${TBB_INCLUDE_DIRS}")
			message(STATUS "TBB_FOUND :        ${TBB_FOUND}")
			message(STATUS "TBB_INCLUDE_DIRS : ${TBB_INCLUDE_DIRS}")
			message(STATUS "TBB_VERSION :      ${TBB_VERSION}")
			message(STATUS "TBB_LIBRARIES :    ${TBB_LIBRARIES}")
		else(TBB_FOUND)
			include(${CMAKE_SOURCE_DIR}/cmake/BuildTBB.cmake)
			set(INSTALL_OWN_TBB 1)
		endif(TBB_FOUND)
	endif(FORCE_OWN_TBB)
endif(ALTCPU)

# This MUST be set before all cuda_add_executable/library calls.
# Otherwise the -Xcompiler argument to NVCC does not contain -fopenmp
# and #pragma omp critical is SILENTLY ignored, giving a wrong binary!!!
# https://github.com/3dem/relion/issues/1016

FIND_PACKAGE( OpenMP REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# -------------------------------------------------------------------------------MPI--
find_package(MPI REQUIRED)

if ("${MPI_CXX_INCLUDE_DIRS}" STREQUAL "")
	if (NOT "${MPI_CXX_INCLUDE_PATH}" STREQUAL "")
		include_directories("${MPI_CXX_INCLUDE_PATH}")
	endif()
else()
	include_directories("${MPI_CXX_INCLUDE_DIRS}")
endif()

message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "MPI_LIBRARIES : ${MPI_LIBRARIES}")
message(STATUS "MPI_CXX_INCLUDE_PATH : ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI_CXX_LIBRARIES : ${MPI_CXX_LIBRARIES}")

message(STATUS "CMAKE_C_COMPILER : ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
message(STATUS "MPI_C_COMPILER : ${MPI_C_COMPILER}")
message(STATUS "MPI_CXX_COMPILER : ${MPI_CXX_COMPILER}")

SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ----------------------------------------------------------Intel Compiler support --
# ----------------------------------  and build flags including MKL and TBB ---------

message(STATUS "CMAKE_CXX_COMPILER_ID : ${CMAKE_CXX_COMPILER_ID}")

if(MKLFFT)
    if (NOT "$ENV{MKLROOT}" STREQUAL "")
        include_directories("$ENV{MKLROOT}/include/fftw")
        message(STATUS "MKL FFTW wrapper header files: $ENV{MKLROOT}/include/fftw")
    else()
        message("COMPILATION MAY FAIL since no MKL FFTW wrapper header files could be found. Please make sure the MKLROOT environmental variable is set.")
    endif()
    add_definitions(-DMKLFFT)
endif(MKLFFT)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel" OR "${CMAKE_CXX_COMPILER}" MATCHES "icpx")
# "Intel" is for classic Intel compiler and "IntelLLVM" is for oneAPI compiler which is supported from CMake 3.20
    if(MKLFFT)
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qopenmp -mkl=parallel -limf ")
        else()	# Intel oneAPI compiler
            SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fiopenmp -qmkl=parallel -limf ")
        endif()
    endif(MKLFFT)
    if(ALTCPU)
        add_definitions(-DFAST_CENTERFFT)
    endif()
else()
    if(MKLFFT)
        # For the time being, let's use the sequential version (as with FFTW)
        link_directories("$ENV{MKLROOT}/lib/intel64")
        SET(FFTW_LIBRARIES mkl_intel_lp64 mkl_sequential mkl_core)
    endif(MKLFFT)
endif()

if(ALTCPU)
    add_definitions(-DALTCPU)
endif()

# ---------------------------------------------------------------USE TEXTURES OR NOT--
if(NOT CudaTexture OR ALTCPU)
    add_definitions(-DPROJECTOR_NO_TEXTURES)
    message(STATUS "Texture interpolation is omitted.")
endif(NOT CudaTexture OR ALTCPU)

# --------------------------------------------------------------------------X11/FLTK--
option(FORCE_OWN_FLTK "FORCE_OWN_FLTK" OFF)
FIND_PACKAGE(X11)
if(GUI)
    if(X11_FOUND)
        set(FLTK_SKIP_FLUID TRUE) #FLUID is not required for relion
        set(FLTK_SKIP_OPENGL TRUE) #OpenGL is not required for relion
        if(NOT FORCE_OWN_FLTK)
            FIND_PACKAGE(FLTK)
		    if(FLTK_FOUND)
				message(STATUS "X11 and FLTK were found")
				message(STATUS "FLTK_LIBRARIES: ${FLTK_LIBRARIES}")
		    else()
				message(STATUS "No FLTK installation was found")
		    endif() 
        endif(NOT FORCE_OWN_FLTK)
        
	if(NOT FLTK_FOUND)
		include(${CMAKE_SOURCE_DIR}/cmake/BuildFLTK.cmake)
		set(INSTALL_OWN_FLTK 1)
	endif(NOT FLTK_FOUND)

    else(X11_FOUND)
        message( STATUS "\n-- ------------------ YOU HAVE NO X11-LIBS ------------------")
        message( STATUS "CCmake found no X11-libs on your system, which are required for the GUI.")
        message( STATUS " You CAN add the flag -DGUI=OFF to avoid using X11" )
        message(FATAL_ERROR "X11 is required for GUI.")
    endif(X11_FOUND)
    
endif(GUI)

# -------------------------------------------------------------------------------FFT--

if(NOT MKLFFT)
	option(FORCE_OWN_FFTW "FORCE_OWN_FFTW" OFF)
	option(AMDFFTW "Use AMD optimized version of FFTW. This needs a new version of GCC (>= 8.3 recommended)." OFF)
	
	set(FFTW_DOUBLE_REQUIRED TRUE)
	set(FFTW_SINGLE_REQUIRED TRUE)
	
	if(AMDFFTW)
		set(FORCE_OWN_FFTW ON)
	endif()

	if(NOT FORCE_OWN_FFTW)
		FIND_PACKAGE(FFTW COMPONENTS SINGLE DOUBLE)
	endif(NOT FORCE_OWN_FFTW)

	if(NOT FFTW_FOUND)
		include(${CMAKE_SOURCE_DIR}/cmake/BuildFFTW.cmake)
	endif(NOT FFTW_FOUND)
endif(NOT MKLFFT)

# ---------------------------------------------------------------------------SIN/COS--

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(sincos    math.h   HAVE_SINCOS)
check_cxx_symbol_exists(__sincos  math.h   HAVE___SINCOS)

if(HAVE_SINCOS)
    add_definitions(-DHAVE_SINCOS)
endif()
if(HAVE___SINCOS)
    add_definitions(-DHAVE___SINCOS)
endif()

# ------------------------------------------------------------------------------TIFF--

find_package(TIFF REQUIRED)
if(TIFF_FOUND)
	add_definitions(-DHAVE_TIFF)
endif()

find_package(ZLIB)
find_package(PNG)
if(PNG_FOUND)
	add_definitions(-DHAVE_PNG)
endif()


# -----------------------------------------------------------------------------Torch--

#option(TORCH "Enable support for LibTorch" OFF)
#if(TORCH)
#    message(STATUS "Torch support requested by user.")
#    include(${CMAKE_SOURCE_DIR}/cmake/BuildTorch.cmake)
#    add_definitions(-D_TORCH_ENABLED)
#endif(TORCH)

option(FETCH_TORCH_MODELS "Do download pre-trained torch models" ON)
if (FETCH_TORCH_MODELS)
    include(${CMAKE_SOURCE_DIR}/cmake/FetchTorchModels.cmake)
endif()

# ----------------------------------------------------------------------COPY SCRIPTS--

if(FORCE_OWN_FFTW)
    install(DIRECTORY external/fftw/lib/ DESTINATION lib FILES_MATCHING PATTERN "*")
endif()

list(APPEND RELION_SCRIPT_FILES
        star_printtable
        star_plottable
        star_loopheader
        star_datablock_stack
        star_datablock_singlefiles
        star_datablock_ctfdat
        qsub.csh
        it.py
        schemegui.py
        class_ranker.py
)

add_custom_target(copy_scripts ALL)
                                
foreach (SCRIPT_FILE ${RELION_SCRIPT_FILES})
    add_custom_command(TARGET copy_scripts POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E
                     copy ${CMAKE_SOURCE_DIR}/scripts/${SCRIPT_FILE} 
                     ${CMAKE_BINARY_DIR}/bin/relion_${SCRIPT_FILE} )
endforeach()

install( DIRECTORY ${CMAKE_BINARY_DIR}/bin
	 DESTINATION ${CMAKE_INSTALL_PREFIX}
         USE_SOURCE_PERMISSIONS
         FILES_MATCHING PATTERN "*")

# install fltk if we built our own version
if(INSTALL_OWN_FLTK)
    install(DIRECTORY external/fltk/lib/ DESTINATION lib FILES_MATCHING PATTERN "*")
endif()

# -----------------------------------------------------------------RELION COMPONENTS--
option(BUILD_SHARED_LIBS "BUILD_SHARED_LIBS" OFF)
message("BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}")
if(BUILD_SHARED_LIBS)
    message(STATUS "Building shared libs (smaller build size and binaries)")
else()
    message(STATUS "Building static libs (larger build size and binaries)")    
endif()


ADD_SUBDIRECTORY(src/apps)

#message(STATUS "CUDA option = ${CUDA}")
#message(STATUS "ALTCPU option = ${ALTCPU}")
#message(STATUS "DoublePrec_CPU option = ${DoublePrec_CPU}")
#message(STATUS "DoublePrec_ACC option = ${DoublePrec_ACC}")
#message(STATUS "MKLFFT option = ${MKLFFT}")
#message(STATUS "CudaTexture option = ${CudaTexture}")

#get_directory_property( DirDefs COMPILE_DEFINITIONS )
#message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}" )
#message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")
#message(STATUS "CMAKE_C_COMPILER : ${CMAKE_C_COMPILER}")
#message(STATUS "CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
#message(STATUS "MPI_C_COMPILER : ${MPI_C_COMPILER}")
#message(STATUS "MPI_CXX_COMPILER : ${MPI_CXX_COMPILER}")
#message(STATUS "CMAKE_EXE_LINKER_FLAGS : ${CMAKE_EXE_LINKER_FLAGS}")

# -----------------------------------------------------------------------------TESTS--
# Include testing flag(s) as precomiler 
# definitions and include test directives
#enable_testing()
#include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)

option(BUILD_TESTS "Build and configure tests" OFF)
if(BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
endif()

# ----------------------------------------------------------PRINT OUT ALL CMAKE VARS--
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

