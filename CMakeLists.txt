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

# ----------------------------------------------------------------SET CXX STANDARD--
set(CXX_STANDARD_REQUIRED ON)
if(SYCL)
	set(CMAKE_CXX_STANDARD 17)
else()
	set(CMAKE_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_EXTENSIONS OFF)

# ------------------OPTIONS WHICH ARE NEEDED TO SET BUILD-TYPES (COMPILATION FLAGS)--
# -------------------------------------------------------------------FURTHER OPTIONS--

# CUDA on by default, so check for other GPU/CPU-accelration request and possible conflicting dual-requests
option(HIP "Enable HIP GPU acceleration" OFF)
option(SYCL "Enable SYCL acceleration" OFF)
option(ALTCPU "Enable Accelerated CPU version" OFF)
if(HIP)
	set(HIP_ARCH "" CACHE STRING "HIP architecture to use")
	if (HIP_ARCH STREQUAL "")
		message(STATUS "Setting fallback HIP_ARCH=gfx906")
		set(AMDGPU_TARGETS "gfx906" CACHE STRING
				"Semicolon delimited list of AMD GPU targets to compile for (e.g gf908;gfx90a)")
	else ()
		message(STATUS "Using provided HIP_ARCH=${HIP_ARCH}")
		set(AMDGPU_TARGETS "${HIP_ARCH}" CACHE STRING
				"Semicolon delimited list of AMD GPU targets to compile for (e.g gf908;gfx90a)")
	endif ()
	option(CUDA "Enable CUDA GPU acceleration" OFF)
	if(HIP AND CUDA)
		message(FATAL_ERROR "You cannot build with both CUDA=ON and HIP=ON.  Please choose one and rerun CMAKE")
	endif()
elseif(SYCL)
  option(CUDA "Enable CUDA GPU acceleration" OFF)
  if(SYCL AND CUDA)
    message(FATAL_ERROR "You cannot build with both CUDA=ON and SYCL=ON. Please choose one and rerun CMAKE")
  endif()

  if(DEFINED SYCL_CUDA_COMPILE)
	if(DEFINED SYCL_CUDA_TARGET)
	  message(STATUS "Using provided SYCL_CUDA_TARGET: sm_${SYCL_CUDA_TARGET}")
	else(DEFINED SYCL_CUDA_TARGET)
	  message(STATUS "Setting fallback SYCL_CUDA_TARGET: sm_50")
	  set(SYCL_CUDA_TARGET "50")
	endif(DEFINED SYCL_CUDA_TARGET)
  endif(DEFINED SYCL_CUDA_COMPILE)
  if(DEFINED SYCL_HIP_COMPILE)
	if(DEFINED SYCL_HIP_TARGET)
	  message(STATUS "Using provided SYCL_HIP_TARGET: ${SYCL_HIP_TARGET}")
	else(DEFINED SYCL_HIP_TARGET)
	  message(STATUS "Setting fallback SYCL_HIP_TARGET: gfx906")
	  set(SYCL_HIP_TARGET "gfx906")
	endif(DEFINED SYCL_HIP_TARGET)
  endif(DEFINED SYCL_HIP_COMPILE)
elseif(ALTCPU)
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
option(DeviceTexture "Texture Memory on the Device" ON)

if(ALTCPU)
    message(STATUS "ALTCPU enabled - Building CPU-accelerated version of RELION")
elseif(SYCL)
    message(STATUS "SYCL enabled - Building SYCL-accelerated version of RELION")
elseif(HIP)
    message(STATUS "HIP enabled - Building HIP-accelerated version of RELION")
elseif(CUDA)
    message(STATUS "CUDA enabled - Building CUDA-accelerated version of RELION")
else()
    message(STATUS "Acceleration not enabled - Building CPU-only version of RELION")
endif()

if(CUDA OR ALTCPU)
    add_definitions(-DACC_HIP=3 -DACC_CUDA=2 -DACC_CPU=1)
endif()

# -----------------------------------------------DOUBLE PRECISION (GPU-CODE) OR NOT--
if(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to double")
else(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to single")
    add_definitions(-DRELION_SINGLE_PRECISION)
endif(DoublePrec_CPU)

if(DoublePrec_ACC)
    message(STATUS "Setting accelerated code precision to double")
    add_definitions(-DACC_DOUBLE_PRECISION)
	  set(DeviceTexture FALSE)
else(DoublePrec_ACC)
    message(STATUS "Setting accelerated code precision to single")
endif(DoublePrec_ACC)

if(MDT_TYPE_CHECK)
    message("-- MetaDataTable type check is enabled.")
    add_definitions(-DMETADATA_TABLE_TYPE_CHECK)
endif()

if(CUDA)
    # -----------------------------------------------------------------------------CUDA--
    # DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html
    FIND_PACKAGE(CUDA)

    # "-arch=sm_XX" alone is equal to "-arch=compute_XX -code=compute_XX,sm_XX"
    # compute_XX generates PTX; thus GPUs newer than XX can also use the binary.
    # See https://docs.nvidia.com/cuda/cuda-compiler-driver-nvcc/#gpu-architecture-arch.
    set(CUDA_ARCH "" CACHE STRING "CUDA architecture to use")
    if (CUDA_ARCH STREQUAL "")
        if (${CUDA_VERSION} VERSION_LESS "13.0") # Pascal (e.g. 1080)
            message(STATUS "Setting fallback CUDA_ARCH=50 (Pascal)")
            set(CUDARCH "-arch=sm_50")
        else() # CUDA >= 13.0 no longer supports Pascal. 75 = Turing (e.g. 2080) 
            message(STATUS "Setting fallback CUDA_ARCH=75 (Turing)")
            set(CUDARCH "-arch=sm_75")
        endif()
    else ()
        message(STATUS "Using provided CUDA_ARCH=${CUDA_ARCH}")
        set(CUDARCH "-arch=sm_${CUDA_ARCH}")
    endif ()

    if(CUDA_FOUND)
      message(STATUS "Using cuda wrapper to compile....")
      if( (NOT ${CUDA_VERSION} VERSION_LESS "7.5") AND (NOT DoublePrec_ACC) )
          message(STATUS "Cuda version is >= 7.5 and single-precision build, enable double usage warning.")
		      set(WARN_DBL "--ptxas-options=-warn-double-usage") # cuda>=7.5
	    elseif( ${CUDA_VERSION} VERSION_LESS "7.0")
		      message(WARNING "Cuda version is less than 7.0, so relion will be compiled without GPU support.")
          message(STATUS "Using non-cuda compilation....")
		      set(CUDA OFF)
      endif()
      if(CUDA)
          add_definitions(-D_CUDA_ENABLED)
      endif()
    else(CUDA_FOUND)
        message(FATAL_ERROR "CUDA enabled but unlable to locate packages...")
    endif(CUDA_FOUND)
elseif(HIP)
    # ------------------------------------------------------------------------------HIP--
    if (DEFINED ENV{ROCM_PATH})
        include($ENV{ROCM_PATH}/lib/cmake/hip/FindHIP.cmake)
    elseif(DEFINED ENV{HIP_PATH})
        include($ENV{HIP_PATH}/cmake/FindHIP.cmake)
    elseif(EXISTS "${CMAKE_SOURCE_DIR}/cmake/FindHIP.cmake")
        include(${CMAKE_SOURCE_DIR}/cmake/FindHIP.cmake)
    else()
        include(/opt/rocm/hip/cmake/FindHIP.cmake)
    endif()

    find_package(HIP REQUIRED)

    if(HIP_FOUND)
      message(STATUS "Using hip wrapper to compile....")
      if( ${HIP_VERSION} VERSION_LESS "5.0")
		      message(WARNING "ROCm version is less than 5.0, so relion will be compiled without GPU support.")
          set(HIP OFF)
      endif()
      if(HIP)
          add_definitions(-D_HIP_ENABLED)
          foreach(hiplib hipfft hipcub rocfft rocprim hiprand rocrand)
            find_package(${hiplib} REQUIRED)
            if(${hiplib}_FOUND)
                message(STATUS "Found ${${hiplib}_LIBRARIES} in ${${hiplib}_INCLUDE_DIR}")
            endif()
          endforeach()
      endif()
    else(HIP_FOUND)
        message(FATAL_ERROR "HIP enabled but unlable to locate packages. ROCm >= 5.0 is required to configure RELION with HIP using CMake.")
    endif(HIP_FOUND)
elseif(SYCL)
	add_definitions(-D_SYCL_ENABLED=1 -DUSE_MPI_COLLECTIVE)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER}" MATCHES "icpx")
		add_definitions(-D_DPCPP_ENABLED=1)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fhonor-infinities -fhonor-nans -fsigned-zeros")
	endif()
elseif(ALTCPU)
	add_definitions(-DALTCPU=1 -DUSE_MPI_COLLECTIVE)
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fhonor-infinities -fhonor-nans -fsigned-zeros")
	endif()
endif()
 
# ----------------------------------------------------------INCLUDE ALL BUILD TYPES---
# This has to be AFTER project() and find_package(HIP).
# Up to RELION 5.0.0, this must come BEFORE find_package(CUDA) as discussed in issue #1032.
# Takanori fixed this in 5.0.1.
include(${CMAKE_SOURCE_DIR}/cmake/BuildTypes.cmake)

# ------------------------------------------------------------------ALLOCATOR CHOICE--
option(CachedAlloc "CachedAlloc" ON)
if(NOT CachedAlloc)
    if (CUDA)
        add_definitions(-DCUDA_NO_CUSTOM_ALLOCATION)
    endif()
    if (HIP)
        add_definitions(-DHIP_NO_CUSTOM_ALLOCATION)
    endif()
	message(STATUS "Cached allocation is disabled.")
endif(NOT CachedAlloc)
option(CustomAllocMemGuards "CustomAllocMemGuards" OFF)
if(CustomAllocMemGuards)
    add_definitions(-DCUSTOM_ALLOCATOR_MEMGUARD)
    message(STATUS "Abort on out of bound write.")
endif(CustomAllocMemGuards)
# -------------------------------------------------------------FORCE USE OF STL-LIBS--
if (CUDA)
    option(CudaForceSTL "CudaForceSTL" OFF)
    if(CudaForceSTL)
        add_definitions(-DCUDA_FORCESTL)
        message(STATUS "Building cuda files using stl-libs for sort, min and max.")
    endif(CudaForceSTL)
elseif (HIP)
    option(HipForceSTL "HipForceSTL" OFF)
    if(HipForceSTL)
        add_definitions(-DHIP_FORCESTL)
        message(STATUS "Building hip files using stl-libs for sort, min and max.")
    endif(HipForceSTL)
elseif (SYCL)
	option(SyclForceOneDPL "SyclForceOneDPL" OFF)
	if(SyclForceOneDPL)
		add_definitions(-DUSE_ONEDPL)
		if(EXISTS "$ENV{DPL_ROOT}/linux/include/oneapi/dpl/algorithm")
			include_directories("$ENV{DPL_ROOT}/linux/include")
			message(STATUS "Building SYCL files using oneDPL: $ENV{DPL_ROOT}/linux/include")
		else()
			message(WARNING "oneDPL header is not found. Please set DPL_ROOT environment or include oneDPL header directory for CMAKE_CXX_FLAGS")
		endif()
	endif(SyclForceOneDPL)
endif()
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
if(HIP)
    set(CMAKE_CXX_COMPILER hipcc)
    message(STATUS "CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
elseif(SYCL)
#	set(CMAKE_C_COMPILER icx)
#	set(CMAKE_CXX_COMPILER icpx)
	message(STATUS "CMAKE_C_COMPILER : ${CMAKE_C_COMPILER}")
	message(STATUS "CMAKE_CXX_COMPILER : ${CMAKE_CXX_COMPILER}")
endif()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ----------------------------------------------------------------- OpenMP ----------
# This MUST be set before all cuda_add_executable/library calls.
# Otherwise the -Xcompiler argument to NVCC does not contain -fopenmp
# and #pragma omp critical is SILENTLY ignored, giving a wrong binary!!!
# https://github.com/3dem/relion/issues/1016

if(SYCL)    # Intel DPC++ specific
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fiopenmp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fiopenmp")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fiopenmp")

	set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-fPIC -std=c++17 ${CMAKE_CXX_FLAGS}")
else()
	FIND_PACKAGE(OpenMP REQUIRED)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

	set(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "-fPIC -std=c++17 ${CMAKE_CXX_FLAGS}")
endif()

# ----------------------------------------------------------Intel Compiler support --
# ----------------------------------  and build flags including MKL and IPP ---------

message(STATUS "CMAKE_CXX_COMPILER_ID : ${CMAKE_CXX_COMPILER_ID}")

if(MKLFFT)
	if (NOT "$ENV{MKLROOT}" STREQUAL "")
		include_directories("$ENV{MKLROOT}/include/fftw")
		message(STATUS "MKL FFTW wrapper header files: $ENV{MKLROOT}/include/fftw")
	else()
		message(WARNING "COMPILATION MAY FAIL since no MKL FFTW wrapper header files could be found. Please make sure the MKLROOT environmental variable is set.")
	endif()
	add_definitions(-DMKLFFT)
endif(MKLFFT)

# NVCC calls ICC several times. The first call is to check compiler compatibility
# and -Xcompiler is not passed. Thus, even with -diag-disable=10441,10412,
# the deprecation warning cannot be completely silenced.
# https://github.com/3dem/relion/issues/1038#issuecomment-1831312535
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel" OR "${CMAKE_CXX_COMPILER}" MATCHES "icpx")
# "Intel" is for classic Intel compiler and "IntelLLVM" is for oneAPI compiler which is supported from CMake 3.20
	if(ALTCPU OR SYCL)
		add_definitions(-DFAST_CENTERFFT)
	endif()

	if(MKLFFT)
		if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qopenmp -qmkl=parallel -limf")
		else()
			SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -diag-disable=10441,10412")
			SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qopenmp -mkl=parallel -limf -diag-disable=10441,10412")
		endif()
	endif(MKLFFT)

	if(INTEL_IPP OR DEFINED ENV{IPPROOT})
		if(EXISTS "$ENV{IPPROOT}/include/ipp.h" AND (EXISTS "$ENV{IPPROOT}/lib/libipps.a" OR EXISTS "$ENV{IPPROOT}/lib/intel64/libipps.a"))
			SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I$ENV{IPPROOT}/include")
			if(EXISTS "$ENV{IPPROOT}/lib/libipps.a")
				SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L$ENV{IPPROOT}/lib -lipps")
			else()
				SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L$ENV{IPPROOT}/lib/intel64 -lipps")
			endif()
#			if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
#				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qipp")
#				SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -qipp -qipp-link=static")
#			else()
#				SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -diag-disable=10441,10412")
#				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipp -diag-disable=10441,10412")
#				SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ipp -ipp-link=static")
#			endif()
			add_definitions(-DUSE_IPP)
			message(STATUS "Using Intel IPP library: $ENV{IPPROOT}")
		endif()
	endif()
else()
	if(MKLFFT)
		# For the time being, let's use the sequential version (as with FFTW)
		link_directories("$ENV{MKLROOT}/lib/intel64")
		SET(FFTW_LIBRARIES mkl_intel_lp64 mkl_sequential mkl_core)
	endif(MKLFFT)
endif()

# ---------------------------------------------------------------USE TEXTURES OR NOT--
if(NOT DeviceTexture OR ALTCPU OR SYCL)
    add_definitions(-DPROJECTOR_NO_TEXTURES)
    message(STATUS "Texture interpolation is omitted.")
endif()

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

# ------------------------------------------------------------------ZLIB, PNG & JPEG--
find_package(ZLIB)
find_package(PNG)
if(PNG_FOUND)
	add_definitions(-DHAVE_PNG)
endif(PNG_FOUND)

find_package(JPEG)
if(JPEG_FOUND)
	add_definitions(-DHAVE_JPEG)
endif(JPEG_FOUND)

# -----------------------------------------------------------------PYTHON DEPENDENCIES--


set(CONDA_ENV_NAME "relion-5.0" CACHE STRING "The Conda environment name")
include(${CMAKE_SOURCE_DIR}/cmake/FindCondaPython.cmake)

list(APPEND RELION_PYTHON_WRAPPERS
        python_classranker
        python_blush
        python_dynamight
        python_topaz
        python_modelangelo
        python_fetch_weights
	python_tomo_import
	python_tomo_exclude_tilt_images
	python_tomo_align_tilt_series
	python_tomo_view
	python_tomo_pick
	python_tomo_get_particle_poses
	python_tomo_denoise
	filament_selection
)

foreach (SCRIPT_FILE ${RELION_PYTHON_WRAPPERS})
    configure_file(
            ${CMAKE_SOURCE_DIR}/scripts/${SCRIPT_FILE}.in
            ${CMAKE_BINARY_DIR}/bin/relion_${SCRIPT_FILE}
            @ONLY
    )
endforeach()

option(FETCH_WEIGHTS "Download model weights for dependent packages" ON)
if (FETCH_WEIGHTS)
	message(STATUS
			"Will attempts to download model weights for dependent packages..."
			"\n"
			"NOTE: You can skip this step by setting -DFETCH_WEIGHTS=OFF"
	)

	execute_process(
			COMMAND bash ${CMAKE_BINARY_DIR}/bin/relion_python_fetch_weights
			WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
			RESULT_VARIABLE SCRIPT_RESULT
	)

	# Check the result of the script execution
	if(SCRIPT_RESULT EQUAL 0)
		message(STATUS "Successfully downloaded model weights")
	else()
		message(WARNING
				"Failed to download model weights. YOU CAN IGNORE THIS. "
				"The downloads will be attempted again when the model weights are required later. "
				"\n"
				"NOTE: You can skip this step by setting -DFETCH_WEIGHTS=OFF"
		)
	endif()
else()
	message(STATUS "Omitting download of model weights for dependent packages")
endif ()

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
#message(STATUS "BUILD_SHARED_LIBS = ${BUILD_SHARED_LIBS}")
if(BUILD_SHARED_LIBS)
    message(STATUS "Building shared libs (smaller build size and binaries)")
else()
    message(STATUS "Building static libs (larger build size and binaries)")
endif()


ADD_SUBDIRECTORY(src/apps)

#message(STATUS "CUDA option = ${CUDA}")
#message(STATUS "ALTCPU option = ${ALTCPU}")
#message(STATUS "HIP option = ${HIP}")
#message(STATUS "SYCL option = ${SYCL}")
#message(STATUS "DoublePrec_CPU option = ${DoublePrec_CPU}")
#message(STATUS "DoublePrec_ACC option = ${DoublePrec_ACC}")
#message(STATUS "MKLFFT option = ${MKLFFT}")
#message(STATUS "DeviceTexture option = ${DeviceTexture}")

#get_directory_property( DirDefs COMPILE_DEFINITIONS )
#message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}" )
message(STATUS "CMAKE_C_FLAGS : ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS : ${CMAKE_CXX_FLAGS}")
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
