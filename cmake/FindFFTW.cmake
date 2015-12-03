# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#

set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")

if(DoublePrec_CPU)
   # set fftw lib to use double precision
	set(ext_conf_flags_fft --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	set(libfft "fftw3")
else(DoublePrec_CPU)
	# set fftw lib to use single precision
	set(ext_conf_flags_fft --enable-shared --enable-float --prefix=${FFTW_EXTERNAL_PATH})
	set(libfft "fftw3f")
endif(DoublePrec_CPU)	

## ------------------------------------------------------------------- SYSTEM LIBS? --
if(NOT OWN_FFTW)
    set(LIB_PATHFFT $ENV{FFTW_LIB})
    set(INC_PATHFFT $ENV{FFTW_INCLUDE})
 
    find_library(FFTW_LIBRARIES  NAMES ${libfft}  PATHS ${LIB_PATHFFT} )  

    message(STATUS "Looking for fft header ...")
    if(DEFINED ENV{FFTW_INCLUDE})
        find_path(FFTW_PATH     NAMES fftw3.h  PATHS ${INC_PATHFFT} )
        find_path(FFTW_INCLUDES NAMES fftw3.h  PATHS ${INC_PATHFFT} )
    else()
        find_path(FFTW_PATH         NAMES fftw3.h )
        find_path(FFTW_INCLUDES     NAMES fftw3.h )
    endif()
        
    find_library(FFTW_LIBRARIES PATHS $ENV{FFTW_LIB} $ENV{FFTW_HOME} )
    
    if(FFTW_PATH)
        set(FFTW_FOUND TRUE)
        message( STATUS "found system-wide installed fft lib") 
    endif()
    
    if(FFTW_FIND_REQUIRED AND (NOT FFTW_FOUND))
        message( STATUS "\n-- ------------------ YOU HAVE NO FFTW-LIBS ------------------")
        message( STATUS "CCmake found no fftw-libs on your system. Either ")
        message( STATUS "     a) Make sure cmake can find them ")
        message( STATUS "     b) Add the flag -DOWN_FFT=ON to your cmake command to build a local fftw-lib")
        message( STATUS " We recommend making use of your system lib")
        message( STATUS "--------------------------------------------------------")
    	message( FATAL_ERROR "FFTW is required." )	
    endif(FFTW_FIND_REQUIRED AND (NOT FFTW_FOUND))

else()  ## ---------------------------------------------------------------- BUILD EXT LIBS? --

   include(${CMAKE_SOURCE_DIR}/cmake/BuildFFTW.cmake)
        
endif() 
message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")