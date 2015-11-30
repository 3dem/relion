# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#
#
# Relion (to be of any expedient use) neess MPI (thread) 
# support, so we will _require_ thread-enabled fftw.
#
# Double precision is default, single precision can be 
# opted into by specifying this in CMakeLists.txt

# PRECISION OPTION
if(DoublePrec_CPU)
	# set fftw lib to use single (f=float) precision
	set(fft "fftw3")
else(DoublePrec_CPU)
	# set fftw lib to use double precision
	set(fft "fftw3f")
endif(DoublePrec_CPU)	

set(LIB_PATHFFT $ENV{FFTW_LIB})
set(INC_PATHFFT $ENV{FFTW_INCLUDE})

find_library(FFTW_LIBRARIES  NAMES ${fft}  PATHS ${LIB_PATHFFT})

if(DEFINED ENV{FFTW_INCLUDE})
    find_path(FFTW_PATH     NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
    find_path(FFTW_INCLUDES NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
else()
    find_path(FFTW_PATH         NAMES fftw3.h )
    find_path(FFTW_INCLUDES     NAMES fftw3.h )
endif()

message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")

if(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)
	message(STATUS "${fft} found.")
   set(FFTW_FOUND TRUE)
else(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)
	message(STATUS "${fft} not found.")
   set(FFTW_FOUND FALSE)
endif(FFTW_PATH AND FFTW_INCLUDES AND FFTW_LIBRARIES)

if(FFTW_FIND_REQUIRED AND NOT FFTW_FOUND)
	message( FATAL_ERROR "FFTW package is required." )
endif(FFTW_FIND_REQUIRED AND NOT FFTW_FOUND)
