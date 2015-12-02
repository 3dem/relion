# Double precision is default, single precision can be 
# opted into by specifying this in CMakeLists.txt

set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")
# PRECISION OPTION SPECIFIERS FOR EXTERNAL BUILD.     
 
if(DoublePrec_CPU)
   # set fftw lib to use double precision
	set(ext_conf_flags_fft --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	message(STATUS "ext_conf_flags_fft: ${ext_conf_flags_fft}")
	set(libfft "fftw3")
else(DoublePrec_CPU)
	# set fftw lib to use single precision
	set(ext_conf_flags_fft --enable-shared --enable-float --prefix=${FFTW_EXTERNAL_PATH})
	set(libfft "fftw3f")
endif(DoublePrec_CPU)	

## ------------------------------------------------------------- PREVIOUS EXT LIBS? --

find_path(FFTW_PATH         NAMES fftw3.h         PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH)
find_path(FFTW_INCLUDES     NAMES fftw3.h         PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
find_library(FFTW_LIBRARIES NAMES lib${libfft}.so PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)   

if(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)
    set(FFTW_FOUND TRUE)
    message( STATUS "found previously built external (non-system) fft lib")
endif()
## ----------------------------------------------------------------- NEW EXT LIBS? --  
 
if(NOT FFTW_FOUND)

    include(ExternalProject)
    set(FFTW_EXTERNAL_LIBS_TAR_DIRECTORY  ${FFTW_EXTERNAL_PATH})
    set(FFTW_EXTERNAL_LIBS_EXTRACT_TARGET ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY})
    message( STATUS "no fftw found, the following paths are set for libs/headers TO BE built")
    set(FFTW_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
    set(FFTW_FFTW3_LIB_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
    set(FFTW_FFTW3_BUILD_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)
    
    set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external")
   
    set(CMAKE_INSTALL_PREFIX  ${FFTW_EXTERNAL_PATH})
    externalproject_add(FFTW3
    URL ${FFTW_FFTW3_TAR_FILE}
    URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
    DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
    SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}
    INSTALL_DIR ${FFTW_EXTERNAL_PATH}/fftw3
    BUILD_COMMAND ${MAKE}
#	LOG_CONFIGURE
#	LOG_BUILD
    LOG_INSTALL)

    set(NEW_OWN_FFTW TRUE)
    set(FFTW_FOUND TRUE)
    
    set(FFTW_LIBRARIES "${FFTW_EXTERNAL_PATH}/lib/lib${libfft}.so" )  
    set(FFTW_PATH      "${FFTW_EXTERNAL_PATH}/includes/fftw3.h" )  
    set(FFTW_INCLUDES  "${FFTW_EXTERNAL_PATH}/includes/fftw3.h" )

endif()

