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
if(SINGLE_RELION)
	# set fftw lib to use single (f=float) precision
	set(fftw "fftw3f")
else(SINGLE_RELION)
	# set fftw lib to use double precision
	set(fftw "fftw3")
endif(SINGLE_RELION)	

find_library(FFTW_LIBRARIES  NAMES "${fftw}"  PATHS $ENV{FFTW_LIB})

# PARALLELISM OPTIONS
if(NOT NOTHREAD_RELION)
    find_library(FFTW_THREAD_LIBS  NAMES "fftw3_threads"  PATHS $ENV{FFTW_LIB})	
	list(APPEND FFTW_LIBRARIES ${FFTW_THREAD_LIBS} )
endif(NOT NOTHREAD_RELION)

message(STATUS "Looking for ${fftw}.h ...")
if(DEFINED ENV{FFTW_INCLUDE})
    find_path(FFTW_PATH     NAMES ${fftw}.h  PATHS $ENV{FFTW_INCLUDE} NO_DEFAULT_PATH)
    find_path(FFTW_INCLUDES NAMES ${fftw}.h  PATHS $ENV{FFTW_INCLUDE} NO_DEFAULT_PATH)
else()
    find_path(FFTW_PATH     NAMES ${fftw}.h )
    find_path(FFTW_INCLUDES     NAMES ${fftw}.h )
endif()

#find_library(FFTW_LIBRARIES /opt/tcbsys/fftw/3.3.4-sse2-avx/lib )

message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
if(FFTW_PATH)
   set(FFTW_FOUND TRUE)
endif(FFTW_PATH)
