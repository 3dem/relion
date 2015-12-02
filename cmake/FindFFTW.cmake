# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#

## ------------------------------------------------------------------- SYSTEM LIBS? --
if(NOT FORCE_OWN_FFTW)
    set(LIB_PATHFFT $ENV{FFTW_LIB})
    set(INC_PATHFFT $ENV{FFTW_INCLUDE})
 
    find_library(FFTW_LIBRARIES  NAMES fftw3f  PATHS ${LIB_PATHFFT})  
    find_library(FFTWD_LIBRARIES   NAMES fftw3   PATHS ${LIB_PATHFFT}) 
    list(APPEND FFTW_LIBRARIES ${FFTWD_LIBRARIES} )

    message(STATUS "Looking for fft header ...")
    if(DEFINED ENV{FFTW_INCLUDE})
        find_path(FFTW_PATH     NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
        find_path(FFTW_INCLUDES NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
    else()
        find_path(FFTW_PATH         NAMES fftw3.h )
        find_path(FFTW_INCLUDES     NAMES fftw3.h )
    endif()
        
    find_library(FFTW_LIBRARIES /opt/tcbsys/fftw/3.3.4-sse2-avx/lib )
    
endif()

## ---------------------------------------------------------------- BUILD EXT LIBS? --
if(FFTW_PATH)
    set(FFTW_FOUND TRUE)
    message( STATUS "found system-wide installed fft lib")
else()
    if(ALLOW_OWN_FFTW)
        include(${CMAKE_SOURCE_DIR}/cmake/BuildFFTW.cmake)
    endif()
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

message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")