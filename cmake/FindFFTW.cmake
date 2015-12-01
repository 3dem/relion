# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#
#
# Relion (to be of any expedient use) neess MPI (thread) 
# support, so we will _require_ thread-enabled fftw.
#
# Double precision is default, single precision can be 
# opted into by specifying this in CMakeLists.txt

set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external")
# PRECISION OPTION SPECIFIERS FOR EXTERNAL BUILD.     
 
if(DoublePrec_CPU)
   # set fftw lib to use double precision
	set(ext_conf_flags --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	message(STATUS "ext_conf_flags: ${ext_conf_flags}")
	set(libfft "fftw3")
else(DoublePrec_CPU)
	# set fftw lib to use single precision
	set(ext_conf_flags --enable-shared --enable-float --prefix=${FFTW_EXTERNAL_PATH})
	set(libfft "fftw3f")
endif(DoublePrec_CPU)	

set(FORCE_OWN_BUILD TRUE)
## ---------------------------------------------------------------------- CUFFT? --
set(USE_CUFFT FALSE)
if(CUFFT)
    set(LIB_PATHFFT $ENV{CUDA_HOME}/lib64)
    set(INC_PATHFFT $ENV{CUDA_HOME}/include)
    find_library(FFTW_LIBRARIES  NAMES cufftw PATHS ${LIB_PATHFFT})
    find_library(FFT_LIBRARIES   NAMES cufft  PATHS ${LIB_PATHFFT})
    list(APPEND FFTW_LIBRARIES ${FFT_LIBRARIES} )
    set(fft "cufft")
    
else(CUFFT)
    if(NOT FORCE_OWN_BUILD)
        set(LIB_PATHFFT $ENV{FFTW_LIB})
        set(INC_PATHFFT $ENV{FFTW_INCLUDE})
     
        find_library(FFTW_LIBRARIES  NAMES fftw3f  PATHS ${LIB_PATHFFT})  
        find_library(FFTWD_LIBRARIES   NAMES fftw3   PATHS ${LIB_PATHFFT}) 
        list(APPEND FFTW_LIBRARIES ${FFTWD_LIBRARIES} )
        # PARALLELISM OPTIONS
    
        if(NOT NOTHREAD_RELION)
            find_library(FFTW_THREAD_LIBS  NAMES "fftw3_threads"  PATHS $ENV{FFTW_LIB})	
       	list(APPEND FFTW_LIBRARIES ${FFTW_THREAD_LIBS} )
        endif(NOT NOTHREAD_RELION)
    endif()
endif(CUFFT)

    
## ------------------------------------------------------------------- SYSTEM LIBS? --
if(NOT FORCE_OWN_BUILD)
    message(STATUS "Looking for fft header ...")
    if(DEFINED ENV{FFTW_INCLUDE})
        find_path(FFTW_PATH     NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
        find_path(FFTW_INCLUDES NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
    else()
        find_path(FFTW_PATH         NAMES fftw3.h )
        find_path(FFTW_INCLUDES     NAMES fftw3.h )
    endif()
        
    find_library(FFTW_LIBRARIES /opt/tcbsys/fftw/3.3.4-sse2-avx/lib )
    
    message(STATUS "FFTW_PATH: ${FFTW_PATH}")
    message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
    message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
endif()
if(FFTW_PATH)
   set(FFTW_FOUND TRUE)
   message( STATUS "found system-wide installed fft lib")
else()
## ------------------------------------------------------------- PREVIOUS EXT LIBS? --

    find_path(FFTW_PATH         NAMES fftw3.h         PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH)
    find_path(FFTW_INCLUDES     NAMES fftw3.h         PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
    find_library(FFTW_LIBRARIES NAMES lib${libfft}.so PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)   
    
    if(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)
        set(FFTW_FOUND TRUE)
        message( STATUS "found previously built external (non-system) fft lib")
    endif()
## ----------------------------------------------------------------- NEW EXT LIBS? --  
 
    if((NOT FFTW_FOUND) AND ALLOW_OWN_BUILD)
        include(ExternalProject)
        set(REL_EXTERNAL_LIBS_TAR_DIRECTORY  ${CMAKE_SOURCE_DIR}/external)
        set(REL_EXTERNAL_LIBS_EXTRACT_TARGET ${REL_EXTERNAL_LIBS_TAR_DIRECTORY})
        message( STATUS "\n-- ------------------ HEY HEY HEY -----------------------")
        message( STATUS "no fftw found, WILL build (because ALLOW_OWN_BUILD=TRUE)")
        message( STATUS "the following paths are set for libs/headers TO BE built")
        message( STATUS "--------------------------------------------------------")
        set(REL_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
        set(REL_FFTW3_LIB_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
        set(REL_FFTW3_BUILD_DIR ${REL_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)
        
        set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external")
       
        set(CMAKE_INSTALL_PREFIX  ${FFTW_EXTERNAL_PATH})
        externalproject_add(FFTW3
        URL ${REL_FFTW3_TAR_FILE}
        URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
        DOWNLOAD_DIR ${REL_EXTERNAL_LIBS_TAR_DIRECTORY}
        SOURCE_DIR ${REL_FFTW3_LIB_DIR}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags}
        INSTALL_DIR ${FFTW_EXTERNAL_PATH}
        BUILD_COMMAND ${MAKE}
    #	LOG_CONFIGURE
    #	LOG_BUILD
        LOG_INSTALL)
        
        set(BUILD_EXTERNAL_FFTW TRUE)
        set(FFTW_FOUND TRUE)
    endif()
    
    if(FFTW_FOUND)
        set(FFTW_LIBRARIES "${FFTW_EXTERNAL_PATH}/lib/lib${libfft}.so" )  
        set(FFTW_PATH      "${FFTW_EXTERNAL_PATH}/includes/fftw3.h" )  
        set(FFTW_INCLUDES  "${FFTW_EXTERNAL_PATH}/includes/fftw3.h" )
    endif()
endif() 

message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
