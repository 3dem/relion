cmake_minimum_required(VERSION 2.8.12)
set(RELION_CMAKE_MINIMUM_REQUIRED_VERSION "2.8.12")

#
#if(POLICY CMP0048)
#  cmake_policy(SET CMP0048 NEW)
#endif()
project(Relion)

# Use new policy for OS X @rpath
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

# Add the path to the additional Find<module>.cmake files 
# which are included with the distributed RLEION-code
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

add_definitions(-DINSTALL_LIBRARY_DIR=${CMAKE_INSTALL_PREFIX}/lib/)
add_definitions(-DSOURCE_DIR=${CMAKE_SOURCE_DIR}/src/)

# message(STATUS "INSTALL_LIBRARY_DIR set to ${CMAKE_INSTALL_PREFIX}/lib/")
# message(STATUS "SOURCE_DIR set to ${CMAKE_SOURCE_DIR}/src/")

# ------------------------------------------------------------------RPATH SETTINGS--
if(NOT APPLE)
    # use, i.e. don't skip the full RPATH for the build tree
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


    # the RPATH to be used when installing, but only if it's not a system directory
    LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
    IF("${isSystemDir}" STREQUAL "-1")
        SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    ENDIF("${isSystemDir}" STREQUAL "-1")
endif(NOT APPLE)

# ---------------------------------------------------------SET SPECIFIC BUILD TYPE--
if(NOT ${CMAKE_BUILD_TYPE} STREQUAL "")
	string( TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER )
	
	if( ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "none" ) AND
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "release" ) AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "debug" ) AND
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "relwithdebinfo" ) AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "profiling" )  AND 
	    ( NOT ${CMAKE_BUILD_TYPE_LOWER} STREQUAL "benchmarking" ) )
	     message( FATAL_ERROR "CMAKE_BUILD_TYPE : '${CMAKE_BUILD_TYPE}' is not a valid build type. "
		"Valid options are: 'None', 'Release', 'Debug', 'RelWithDebInfo', and 'Profiling'." )
	endif()
	
	message(STATUS "BUILD TYPE set to '${CMAKE_BUILD_TYPE}'")
	SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: 'None', 'Release', 'Debug', 'RelWithDebInfo', and 'Profiling'.")
else()
	SET(CMAKE_BUILD_TYPE "Release")
	message(STATUS "BUILD TYPE set to the default type:  '${CMAKE_BUILD_TYPE}'")
	string( TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_LOWER )
endif()

# ------------------OPTIONS WHICH ARE NEEDED TO SET BUILD-TYPES (COMPILATION FLAGS)--
# ------------------------------------------------------------------------CUDA-ARCH--
if(NOT DEFINED CUDA_ARCH)
    message(STATUS "Setting fallback CUDA_ARCH=35")
    set(CUDARCH "-arch=sm_35")
else(NOT DEFINED CUDA_ARCH)
    message(STATUS "Using provided CUDA_ARCH=${CUDA_ARCH}")
    set(CUDARCH "-arch=sm_${CUDA_ARCH}")
endif(NOT DEFINED CUDA_ARCH)

# -------------------------------------------------------------------FURTHER OPTIONS--
option(CUDA "Enable CUDA GPU acceleration" ON)
option(DoublePrec_CPU "DoublePrec_CPU" ON)
option(DoublePrec_GPU "DoublePrec_GPU" OFF)
option(CudaTexture "CudaTexture" ON)

# -----------------------------------------------DOUBLE PRECISION (CUDA-CODE) OR NOT--
if(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to double")    
else(DoublePrec_CPU)
    message(STATUS "Setting cpu precision to single")
    add_definitions(-DRELION_SINGLE_PRECISION)
endif(DoublePrec_CPU)

if(DoublePrec_GPU)
    message(STATUS "Setting gpu precision to double")
    add_definitions(-DCUDA_DOUBLE_PRECISION)
	set(CudaTexture FALSE)
else(DoublePrec_GPU)
    message(STATUS "Setting gpu precision to single")
endif(DoublePrec_GPU)

# ----------------------------------------------------------INCLUDE ALL BUILD TYPES--
 #This *has* to be AFTER project()
include(${CMAKE_SOURCE_DIR}/cmake/BuildTypes.cmake)

if(CUDA)
    # -----------------------------------------------------------------------------CUDA--
    # DOC: http://www.cmake.org/cmake/help/v3.0/module/FindCUDA.html
    FIND_PACKAGE(CUDA)
endif()

if(CUDA_FOUND)
    message(STATUS "Using cuda wrapper to compile....")
    if( (NOT ${CUDA_VERSION} VERSION_LESS "7.5") AND (NOT DoublePrec_GPU) )
        message(STATUS "Cuda version is >= 7.5 and single-precision build, enable double usage warning.")
		set(WARN_DBL "--ptxas-options=-warn-double-usage") # cuda>=7.5
	elseif( ${CUDA_VERSION} VERSION_LESS "7.0")
		message(WARNING "Cuda version is less than 7.0, so relion will be compiled without GPU support.")
		set(CUDA OFF)
    endif()
    
    if(CUDA)
    	add_definitions(-DCUDA)
    endif()
else(CUDA_FOUND)
    message(STATUS "Using non-cuda compilation....")
endif(CUDA_FOUND)

# ---------------------------------------------------------------USE TEXTURES OR NOT--
if(NOT CudaTexture)
    add_definitions(-DCUDA_NO_TEXTURES)
    message(STATUS "Texture interpolation is omitted.")
endif(NOT CudaTexture)
# ------------------------------------------------------------------ALLOCATOR CHOICE--
option(CachedAlloc "CachedAlloc" ON)
if(NOT CachedAlloc)
    add_definitions(-DCUDA_NO_CUSTOM_ALLOCATION)
    message(STATUS "Cashed allocation is disabled.")
endif(NOT CachedAlloc)
option(CustomAllocMemGuards "CustomAllocMemGuards" OFF)
if(CustomAllocMemGuards)
    add_definitions(-DCUSTOM_ALLOCATOR_MEMGUARD)
    message(STATUS "Abort on out of bound write.")
endif(CustomAllocMemGuards)
# -------------------------------------------------------------FORCE USE OF STL-LIBS--
option(CudaForceSTL "CudaForceSTL" OFF)
if(CudaForceSTL)
    add_definitions(-DCUDA_FORCESTL)
    message(STATUS "Building cuda files wusing stl-libs for sort, min and max.")
endif(CudaForceSTL)

# ------------------------------------------------------------------------GUI OR NOT--
# Skip FLTK/X11-dependent binaries or not
option(GUI "GUI" ON)
if(NOT GUI)
    message(STATUS "Omitting GUI targets as per your request")
endif()

# -------------------------------------------------------------------------------MPI--
find_package(MPI REQUIRED)

include_directories("${MPI_INCLUDE_PATH}")
message(STATUS "MPI_INCLUDE_PATH : ${MPI_INCLUDE_PATH}")
message(STATUS "MPI_LIBRARIES : ${MPI_LIBRARIES}")
message(STATUS "MPI_CXX_INCLUDE_PATH : ${MPI_CXX_INCLUDE_PATH}")
message(STATUS "MPI_CXX_LIBRARIES : ${MPI_CXX_LIBRARIES}")

SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(FORCE_OWN_FLTK "FORCE_OWN_FLTK" OFF)
# --------------------------------------------------------------------------X11/FLTK--
FIND_PACKAGE(X11)
if(GUI)
    if(X11_FOUND)
        set(FLTK_SKIP_OPENGL TRUE) #OpenGL is not required for relion
        if(NOT FORCE_OWN_FLTK)
            FIND_PACKAGE(FLTK)
		    if(FLTK_FOUND)
				message(STATUS "X11 and FLTK were found")
				message(STATUS "FLTK_LIBRARIES: ${FLTK_LIBRARIES}")
		    else()
				message(STATUS "No FLTK installation was found")
		    endif() 
        endif(NOT FORCE_OWN_FLTK)
        
	if(NOT FLTK_FOUND)
            include(${CMAKE_SOURCE_DIR}/cmake/BuildFLTK.cmake)
            set(INSTALL_OWN_FLTK 1)
        endif(NOT FLTK_FOUND)

    else(X11_FOUND)
        message( STATUS "\n-- ------------------ YOU HAVE NO X11-LIBS ------------------")
        message( STATUS "CCmake found no X11-libs on your system, which are required for the GUI.")
        message( STATUS " You CAN add the flag -DGUI=OFF to avoid using X11" )
        message(FATAL_ERROR "X11 is required for GUI.")
    endif(X11_FOUND)
    
endif(GUI)

# ------------------------------------------------------------------------------FFTW--
option(FORCE_OWN_FFTW "FORCE_OWN_FFTW" OFF)

if(NOT FORCE_OWN_FFTW)
	FIND_PACKAGE(FFTW)
endif(NOT FORCE_OWN_FFTW)
	
if(NOT FFTW_FOUND)
	message(STATUS "No FFTW installation was found")
	include(${CMAKE_SOURCE_DIR}/cmake/BuildFFTW.cmake)
endif(NOT FFTW_FOUND)

include(CheckCXXSymbolExists)
check_cxx_symbol_exists(sincos    math.h   HAVE_SINCOS)
check_cxx_symbol_exists(__sincos  math.h   HAVE___SINCOS)

if(HAVE_SINCOS)
    add_definitions(-DHAVE_SINCOS)
endif()
if(HAVE___SINCOS)
    add_definitions(-DHAVE___SINCOS)
endif()

# ----------------------------------------------------------------------COPY SCRIPTS--

if(FORCE_OWN_FFTW)
    install(DIRECTORY external/fftw/lib/ DESTINATION lib FILES_MATCHING PATTERN "*")
endif()

list(APPEND RELION_SCRIPT_FILES star_printtable 
                                star_plottable
                                star_loopheader
                                star_datablock_stack
                                star_datablock_singlefiles
                                star_datablock_ctfdat
                                qsub.csh)

add_custom_target(copy_scripts ALL)
                                
foreach (SCRIPT_FILE ${RELION_SCRIPT_FILES})
    add_custom_command(TARGET copy_scripts POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E
                     copy ${CMAKE_SOURCE_DIR}/scripts/${SCRIPT_FILE} 
                     ${CMAKE_BINARY_DIR}/bin/relion_${SCRIPT_FILE} )
endforeach()
install( DIRECTORY ${CMAKE_BINARY_DIR}/bin
	 DESTINATION ${CMAKE_INSTALL_PREFIX}
         USE_SOURCE_PERMISSIONS
         FILES_MATCHING PATTERN "*")

install(FILES ${CMAKE_SOURCE_DIR}/src/gui_background.xpm DESTINATION lib)

# install fltk if we built our own version
if(INSTALL_OWN_FLTK)
    install(DIRECTORY external/fltk/lib/ DESTINATION lib FILES_MATCHING PATTERN "*")
endif()

# -----------------------------------------------------------------RELION COMPONENTS--
option(BUILD_SHARED_LIBS "BUILD_SHARED_LIBS" ON)
if(BUILD_SHARED_LIBS)
    message(STATUS "Building shared libs (smaller build size and binaries)")
else()
    message(STATUS "Building static libs (larger build size and binaries)")    
endif()


ADD_SUBDIRECTORY(src/apps)

# -----------------------------------------------------------------------------TESTS--
# Include testing flag(s) as precomiler 
# definitions and include test directives
enable_testing()
include(${CMAKE_SOURCE_DIR}/tests/RelionTests.cmake)

# ----------------------------------------------------------PRINT OUT ALL CMAKE VARS--
#get_cmake_property(_variableNames VARIABLES)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

