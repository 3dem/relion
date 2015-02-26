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
if(FFTW_SINGLE)
	# set fftw lib to use single (f=float) precision
	set(fftw "fftw3f")
else(FFTW_SINGLE)
	# set fftw lib to use double precision
	set(fftw "fftw3")
endif(FFTW_SINGLE)	

find_library(FFTW_LIBRARIES "${fftw}")

# PARALLELISM OPTIONS
if(NOT NOTHREAD_RELION)
        find_library(FFTW_THREAD_LIBS "${fftw}_threads")	
	list(APPEND FFTW_LIBRARIES ${FFTW_THREAD_LIBS} )
endif(NOT NOTHREAD_RELION)

find_path(FFTW_INCLUDES ${fftw}.h)
find_path(FFTW_PATH ${fftw}.h /usr/lib/* /usr/share/* /usr/include/* )

if(FFTW_PATH)
   set(FFTW_FOUND TRUE)
endif(FFTW_PATH)
