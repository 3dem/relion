# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#

set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")

set(LIB_PATHFFT $ENV{FFTW_LIB})
set(INC_PATHFFT $ENV{FFTW_INCLUDE})

unset(FFTW_PATH CACHE)
unset(FFTW_INCLUDES CACHE)
unset(FFTW_LIBRARIES CACHE)

if(NOT DoublePrec_CPU AND NOT DoublePrec_ACC)
    find_library(FFTW_LIBRARIES  NAMES fftw3f  PATHS ${LIB_PATHFFT} $ENV{FFTW_LIB} $ENV{FFTW_HOME} ) 
elseif(DoublePrec_CPU AND DoublePrec_ACC)
    find_library(FFTW_LIBRARIES  NAMES fftw3   PATHS ${LIB_PATHFFT} $ENV{FFTW_LIB} $ENV{FFTW_HOME} ) 
else() 
    find_library(FFTW_LIBSINGLE  NAMES fftw3f  PATHS ${LIB_PATHFFT} $ENV{FFTW_LIB} $ENV{FFTW_HOME} )
    find_library(FFTW_LIBDOUBLE  NAMES fftw3   PATHS ${LIB_PATHFFT} $ENV{FFTW_LIB} $ENV{FFTW_HOME} )
    set(FFTW_LIBRARIES ${FFTW_LIBSINGLE} ${FFTW_LIBDOUBLE})
endif()	   

if(DEFINED ENV{FFTW_INCLUDE})
    find_path(FFTW_PATH     NAMES fftw3.h PATHS ${INC_PATHFFT} )
    find_path(FFTW_INCLUDES NAMES fftw3.h PATHS ${INC_PATHFFT} )
else()
    find_path(FFTW_PATH     NAMES fftw3.h )
    find_path(FFTW_INCLUDES NAMES fftw3.h )
endif()


if(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)
   set(FFTW_FOUND TRUE)
endif(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)

if (FFTW_FOUND)
	message(STATUS "Found FFTW")
	message(STATUS "FFTW_PATH: ${FFTW_PATH}")
	message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
	message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
else(FFTW_FOUND)
	message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
	if(DoublePrec_CPU)
		message(STATUS "Double-precision FFTW was NOT found")
	else(DoublePrec_CPU)
		message(STATUS "Single-precision FFTW was NOT found")
	endif(DoublePrec_CPU)	
endif(FFTW_FOUND)

if(FFTW_FIND_REQUIRED AND NOT FFTW_FOUND)
	message( FATAL_ERROR "FFTW is required." )
endif(FFTW_FIND_REQUIRED AND NOT FFTW_FOUND)
