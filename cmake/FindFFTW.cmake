# CMake helper to locate the needed libraries and headers
# for compilation of RELION binaries.
#
#
# Relion (to be of any expedient use) neess MPI (thread) 
# support, so we will _require_ thread-enabled fftw.
#
# Double precision is default, single precision can be 
# opted into by specifying this in CMakeLists.txt


# PRECISION CHOICES
if(SINGLE_RELION)
	# set fftw lib to use single (f=float) precision
	set(fftw "fftw3f")
else()
	# set fftw lib to use double precision
	set(fftw "fftw3")
endif()	

if(NOT NOTHREAD_RELION)
	set(thread_opt "_threads")
        find_library(FFTW_LIB1 "${fftw}${thread_opt}")
endif()
find_library(FFTW_LIB2 "${fftw}")
set(FFTW_LIBRARIES ${FFTW_LIB1} ${FFTW_LIB2} )



find_path(FFTW_INCLUDES ${fftw}.h)
find_path(FFTW_PATH ${fftw}.h /usr/lib/* /usr/share/* /usr/include/* )

if(FFTW_PATH)
   set(FFTW_FOUND TRUE)
endif()

