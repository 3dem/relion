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


#set(USE_CUFFT FALSE)
if(CUFFT)
    set(LIB_PATHFFT $ENV{CUDA_HOME}/lib64)
    set(INC_PATHFFT $ENV{CUDA_HOME}/include)
    find_library(FFTW_LIBRARIES  NAMES cufftw PATHS ${LIB_PATHFFT})
    find_library(FFT_LIBRARIES   NAMES cufft  PATHS ${LIB_PATHFFT})
    list(APPEND FFTW_LIBRARIES ${FFT_LIBRARIES} )
    set(fft "cufft")
    
else(CUFFT)
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
   
endif(CUFFT)    

message(STATUS "Looking for ${fft}.h ...")
if(DEFINED ENV{FFTW_INCLUDE})
    find_path(FFTW_PATH     NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
    find_path(FFTW_INCLUDES NAMES fftw3.h  PATHS ${INC_PATHFFT} NO_DEFAULT_PATH)
else()
    find_path(FFTW_PATH         NAMES fftw3.h )
    find_path(FFTW_INCLUDES     NAMES fftw3.h )
endif()
    

#find_library(FFTW_LIBRARIES /opt/tcbsys/fftw/3.3.4-sse2-avx/lib )

message(STATUS "FFTW_PATH: ${FFTW_PATH}")
message(STATUS "FFTW_INCLUDES: ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
if(FFTW_PATH)
   set(FFTW_FOUND TRUE)
endif(FFTW_PATH)
