include_directories("${CMAKE_SOURCE_DIR}")

message("Running apps/CMakeLists.txt...")

if(NOT MKLFFT)
	include_directories("${FFTW_INCLUDES}")
endif(NOT MKLFFT)
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_GUI_SRC "${CMAKE_SOURCE_DIR}/src/manualpicker.cpp" "${CMAKE_SOURCE_DIR}/src/gui_*.cpp" "${CMAKE_SOURCE_DIR}/src/displayer.cpp")

# AUTOMATIC VERSIONING
# Reference: https://stackoverflow.com/questions/1435953/how-can-i-pass-git-sha1-to-compiler-as-definition-using-cmake

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1_FULL)
set(RELION_VERSION_SUFFIX "")
if (GIT_SHA1_FULL)
   string(SUBSTRING ${GIT_SHA1_FULL} 0 6 GIT_SHA1_PREFIX)
   if (NOT GIT_SHA1_PREFIX STREQUAL "GITDIR")
      message(STATUS "Git commit ID: ${GIT_SHA1_FULL}")
      set(RELION_VERSION_SUFFIX "-commit-${GIT_SHA1_PREFIX}")
   endif()
endif()
configure_file("${CMAKE_SOURCE_DIR}/src/macros.cpp.in" "${CMAKE_BINARY_DIR}/macros.cpp" @ONLY)

if (ALTCPU)
	file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_BINARY_DIR}/macros.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.cpp" "${CMAKE_SOURCE_DIR}/src/acc/cpu/cpu_kernels/*.cpp" )
	file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/acc/*.h" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.h" "${CMAKE_SOURCE_DIR}/src/acc/cpu/cpu_kernels/*.h" )
else()
	file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_BINARY_DIR}/macros.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/acc/*.cpp" )
	file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/acc/*.h" ) 
endif(ALTCPU)

# Remove GUI files from relion_lib
foreach(GUI_SRC_FILE ${REL_GUI_SRC})
	list(REMOVE_ITEM REL_SRC "${GUI_SRC_FILE}")
endforeach()

file(GLOB REL_JAZ_MATH_SRC "${CMAKE_SOURCE_DIR}/src/jaz/math/*.cpp")
file(GLOB REL_JAZ_MATH_H "${CMAKE_SOURCE_DIR}/src/jaz/math/*.h")
file(GLOB REL_JAZ_IMG_PROC_SRC "${CMAKE_SOURCE_DIR}/src/jaz/image/*.cpp")
file(GLOB REL_JAZ_IMG_PROC_H "${CMAKE_SOURCE_DIR}/src/jaz/image/*.h")
file(GLOB REL_JAZ_MESH_PROC_SRC "${CMAKE_SOURCE_DIR}/src/jaz/mesh/*.cpp")
file(GLOB REL_JAZ_MESH_PROC_H "${CMAKE_SOURCE_DIR}/src/jaz/mesh/*.h")
file(GLOB REL_JAZ_UTIL_SRC "${CMAKE_SOURCE_DIR}/src/jaz/util/*.cpp")
file(GLOB REL_JAZ_UTIL_H "${CMAKE_SOURCE_DIR}/src/jaz/util/*.h")
file(GLOB REL_JAZ_TOMO_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/*.cpp")
file(GLOB REL_JAZ_TOMO_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/*.h")
file(GLOB REL_JAZ_TOMO_PROJ_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/projection/*.cpp")
file(GLOB REL_JAZ_TOMO_PROJ_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/projection/*.h")
file(GLOB REL_JAZ_TOMO_LATTICE_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/lattice/*.cpp")
file(GLOB REL_JAZ_TOMO_LATTICE_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/lattice/*.h")
file(GLOB REL_JAZ_TOMO_MEMBRANE_SRC "${CMAKE_SOURCE_DIR}/src/jaz/membrane/*.cpp")
file(GLOB REL_JAZ_TOMO_MEMBRANE_H "${CMAKE_SOURCE_DIR}/src/jaz/membrane/*.h")
file(GLOB REL_JAZ_TOMO_FILAMENT_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/filament/*.cpp")
file(GLOB REL_JAZ_TOMO_FILAMENT_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/filament/*.h")
file(GLOB REL_JAZ_TOMO_MOTION_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/motion/*.cpp")
file(GLOB REL_JAZ_TOMO_MOTION_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/motion/*.h")
file(GLOB REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/motion/modular_alignment/*.cpp")
file(GLOB REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/motion/modular_alignment/*.h")
file(GLOB REL_JAZ_TOMO_MANIFOLD_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/manifold/*.cpp")
file(GLOB REL_JAZ_TOMO_MANIFOLD_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/manifold/*.h")
file(GLOB REL_JAZ_OPTIM_SRC "${CMAKE_SOURCE_DIR}/src/jaz/optimization/*.cpp")
file(GLOB REL_JAZ_OPTIM_H "${CMAKE_SOURCE_DIR}/src/jaz/optimization/*.h")
file(GLOB REL_JAZ_OPTICS_SRC "${CMAKE_SOURCE_DIR}/src/jaz/optics/*.cpp")
file(GLOB REL_JAZ_OPTICS_H "${CMAKE_SOURCE_DIR}/src/jaz/optics/*.h")
file(GLOB REL_JAZ_OPTICS_DUAL_CONTRAST_SRC "${CMAKE_SOURCE_DIR}/src/jaz/optics/dual_contrast/*.cpp")
file(GLOB REL_JAZ_OPTICS_DUAL_CONTRAST_H "${CMAKE_SOURCE_DIR}/src/jaz/optics/dual_contrast/*.h")
file(GLOB REL_JAZ_SEGMENTATION_SRC "${CMAKE_SOURCE_DIR}/src/jaz/segmentation/*.cpp")
file(GLOB REL_JAZ_SEGMENTATION_H "${CMAKE_SOURCE_DIR}/src/jaz/segmentation/*.h")
file(GLOB REL_JAZ_ATOMIC_SRC "${CMAKE_SOURCE_DIR}/src/jaz/atomic/*.cpp")
file(GLOB REL_JAZ_ATOMIC_H "${CMAKE_SOURCE_DIR}/src/jaz/atomic/*.h")
file(GLOB REL_JAZ_PROGRAMS_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/programs/*.cpp")
file(GLOB REL_JAZ_PROGRAMS_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/programs/*.h")
file(GLOB REL_JAZ_DYNAMO_SRC "${CMAKE_SOURCE_DIR}/src/jaz/tomography/dynamo/*.cpp")
file(GLOB REL_JAZ_DYNAMO_H "${CMAKE_SOURCE_DIR}/src/jaz/tomography/dynamo/*.h")

file(GLOB REL_JAZ_SPA_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/*.cpp")
file(GLOB REL_JAZ_SPA_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/*.h")
file(GLOB REL_JAZ_SPA_OPT_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/optimization/*.cpp")
file(GLOB REL_JAZ_SPA_OPT_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/optimization/*.h")
file(GLOB REL_JAZ_SPA_CTF_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/ctf/*.cpp")
file(GLOB REL_JAZ_SPA_CTF_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/ctf/*.h")
file(GLOB REL_JAZ_SPA_MOTION_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/motion/*.cpp")
file(GLOB REL_JAZ_SPA_MOTION_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/motion/*.h")
file(GLOB REL_JAZ_SPA_MATH_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/math/*.cpp")
file(GLOB REL_JAZ_SPA_MATH_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/math/*.h")
file(GLOB REL_JAZ_SPA_IMG_PROC_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/img_proc/*.cpp")
file(GLOB REL_JAZ_SPA_IMG_PROC_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/img_proc/*.h")
file(GLOB REL_JAZ_SPA_TOMO_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/tomo/*.cpp")
file(GLOB REL_JAZ_SPA_TOMO_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/tomo/*.h")
file(GLOB REL_JAZ_SPA_IO_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/io/*.cpp")
file(GLOB REL_JAZ_SPA_IO_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/io/*.h")
file(GLOB REL_JAZ_SPA_PROGRAMS_SRC "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/programs/*.cpp")
file(GLOB REL_JAZ_SPA_PROGRAMS_H "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/programs/*.h")

file(GLOB REL_SPHERICAL_HARMONICS_SRC "${CMAKE_SOURCE_DIR}/src/spherical-harmonics/*.cc")
file(GLOB REL_SPHERICAL_HARMONICS_H "${CMAKE_SOURCE_DIR}/src/spherical-harmonics/*.h")
file(GLOB REL_D3x3_SRC "${CMAKE_SOURCE_DIR}/src/jaz/d3x3/*.c")
file(GLOB REL_D3x3_H "${CMAKE_SOURCE_DIR}/src/jaz/d3x3/*.h")
file(GLOB REL_LBFGS_SRC "${CMAKE_SOURCE_DIR}/src/lbfgs/*.c")
file(GLOB REL_LBFGS_H "${CMAKE_SOURCE_DIR}/src/lbfgs/*.h")

file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")
file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")
file(GLOB RELION_JAZ_SPA_TARGETS "${CMAKE_SOURCE_DIR}/src/jaz/single_particle/apps/*.cpp")
file(GLOB RELION_JAZ_TOMO_TARGETS "${CMAKE_SOURCE_DIR}/src/jaz/tomography/apps/*.cpp")

if(RELION_TEST)
	file(GLOB RELION_JAZ_EXPERIMENTAL_TARGETS "${CMAKE_SOURCE_DIR}/src/jaz/scripts/*.cpp")
endif(RELION_TEST)

set(GUI_TARGETS maingui display manualpick)

set(TORCH_TARGETS class_ranker)

#--Remove apps using X11 if no GUI--
if(NOT GUI)
	foreach(TARGET ${GUI_TARGETS})
		list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/${TARGET}.cpp")
	endforeach()
endif(NOT GUI)

#--Remove apps for testing--

set(TEST_TARGETS movie_reconstruct double_reconstruct_openmp cs_fit ctf_nyquist_test free_aberration_plot split_stack defocus_stats double_bfac_fit interpolation_test motion_diff paper_data_synth Zernike_test vis_delocalisation vis_Ewald_weight mpi_tester)
if(NOT RELION_TEST)
	foreach(TARGET ${TEST_TARGETS})
		list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/${TARGET}.cpp")
	endforeach()
endif(NOT RELION_TEST)

# relion_lib is STATIC or SHARED type based on BUILD_SHARED_LIBS=ON/OFF
# relion_lib only contains non-X11 parts
# relion_gui_lib is where the X11 code is placed
if(BUILD_SHARED_LIBS)
	add_library(relion_lib SHARED
		${REL_SRC} ${REL_SRC_H} ${REL_HP}
		${REL_JAZ_MATH_SRC} ${REL_JAZ_MATH_H}
		${REL_JAZ_IMG_PROC_SRC} ${REL_JAZ_IMG_PROC_H}
		${REL_JAZ_MESH_PROC_SRC} ${REL_JAZ_MESH_PROC_H}
		${REL_JAZ_UTIL_SRC} ${REL_JAZ_UTIL_H}
		${REL_JAZ_TOMO_SRC} ${REL_JAZ_TOMO_H}
		${REL_JAZ_TOMO_PROJ_SRC} ${REL_JAZ_TOMO_PROJ_H}
		${REL_JAZ_TOMO_MEMBRANE_SRC} ${REL_JAZ_TOMO_MEMBRANE_H}
		${REL_JAZ_TOMO_FILAMENT_SRC} ${REL_JAZ_TOMO_FILAMENT_H}
		${REL_JAZ_TOMO_MOTION_SRC} ${REL_JAZ_TOMO_MOTION_H}
		${REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_SRC} ${REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_H}
		${REL_JAZ_TOMO_MANIFOLD_SRC} ${REL_JAZ_TOMO_MANIFOLD_H}
		${REL_JAZ_TOMO_LATTICE_SRC} ${REL_JAZ_TOMO_LATTICE_H}
		${REL_JAZ_OPTIM_SRC} ${REL_JAZ_OPTIM_H}
		${REL_JAZ_OPTICS_SRC} ${REL_JAZ_OPTICS_H}
		${REL_JAZ_OPTICS_DUAL_CONTRAST_SRC} ${REL_JAZ_OPTICS_DUAL_CONTRAST_H}
		${REL_JAZ_SEGMENTATION_SRC} ${REL_JAZ_SEGMENTATION_H}
		${REL_JAZ_ATOMIC_SRC} ${REL_JAZ_ATOMIC_H}
		${REL_JAZ_PROGRAMS_SRC} ${REL_JAZ_PROGRAMS_H}
		${REL_JAZ_DYNAMO_SRC} ${REL_JAZ_DYNAMO_H}
		${REL_JAZ_SPA_SRC} ${REL_JAZ_SPA_H}
		${REL_JAZ_SPA_OPT_SRC} ${REL_JAZ_SPA_OPT_H}
		${REL_JAZ_SPA_CTF_SRC} ${REL_JAZ_SPA_CTF_H}
		${REL_JAZ_SPA_MOTION_SRC} ${REL_JAZ_SPA_MOTION_H}
		${REL_JAZ_SPA_MATH_SRC} ${REL_JAZ_SPA_MATH_H}
		${REL_JAZ_SPA_IMG_PROC_SRC} ${REL_JAZ_SPA_IMG_PROC_H}
		${REL_JAZ_SPA_TOMO_SRC} ${REL_JAZ_SPA_TOMO_H}
		${REL_JAZ_SPA_IO_SRC} ${REL_JAZ_SPA_IO_H}
		${REL_JAZ_SPA_PROGRAMS_SRC} ${REL_JAZ_SPA_PROGRAMS_H}
		${REL_D3x3_SRC} ${REL_D3x3_H}
		${REL_LBFGS_SRC} ${REL_LBFGS_H}
		${REL_SPHERICAL_HARMONICS_SRC} ${REL_SPHERICAL_HARMONICS_H})

	install(TARGETS relion_lib LIBRARY DESTINATION lib)
	if(GUI)
		add_library(relion_gui_lib SHARED ${REL_GUI_SRC} ${REL_SRC_H} ${REL_HP})
		install(TARGETS relion_gui_lib LIBRARY DESTINATION lib)
	endif(GUI)
else()
	add_library(relion_lib STATIC
		${REL_SRC} ${REL_SRC_H} ${REL_HP}
		${REL_JAZ_MATH_SRC} ${REL_JAZ_MATH_H}
		${REL_JAZ_IMG_PROC_SRC} ${REL_JAZ_IMG_PROC_H}
		${REL_JAZ_MESH_PROC_SRC} ${REL_JAZ_MESH_PROC_H}
		${REL_JAZ_UTIL_SRC} ${REL_JAZ_UTIL_H}
		${REL_JAZ_TOMO_SRC} ${REL_JAZ_TOMO_H}
		${REL_JAZ_TOMO_PROJ_SRC} ${REL_JAZ_TOMO_PROJ_H}    
		${REL_JAZ_TOMO_MEMBRANE_SRC} ${REL_JAZ_TOMO_MEMBRANE_H}
		${REL_JAZ_TOMO_FILAMENT_SRC} ${REL_JAZ_TOMO_FILAMENT_H}
		${REL_JAZ_TOMO_MOTION_SRC} ${REL_JAZ_TOMO_MOTION_H}
		${REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_SRC} ${REL_JAZ_TOMO_MOTION_MODULAR_ALIGNMENT_H}
		${REL_JAZ_TOMO_MANIFOLD_SRC} ${REL_JAZ_TOMO_MANIFOLD_H}
		${REL_JAZ_TOMO_LATTICE_SRC} ${REL_JAZ_TOMO_LATTICE_H}    
		${REL_JAZ_OPTIM_SRC} ${REL_JAZ_OPTIM_H}
		${REL_JAZ_OPTICS_SRC} ${REL_JAZ_OPTICS_H}
		${REL_JAZ_OPTICS_DUAL_CONTRAST_SRC} ${REL_JAZ_OPTICS_DUAL_CONTRAST_H}
		${REL_JAZ_SEGMENTATION_SRC} ${REL_JAZ_SEGMENTATION_H} 
		${REL_JAZ_ATOMIC_SRC} ${REL_JAZ_ATOMIC_H} 
		${REL_JAZ_PROGRAMS_SRC} ${REL_JAZ_PROGRAMS_H} 
		${REL_JAZ_DYNAMO_SRC} ${REL_JAZ_DYNAMO_H} 
		${REL_JAZ_SPA_SRC} ${REL_JAZ_SPA_H}
		${REL_JAZ_SPA_OPT_SRC} ${REL_JAZ_SPA_OPT_H}
		${REL_JAZ_SPA_CTF_SRC} ${REL_JAZ_SPA_CTF_H}
		${REL_JAZ_SPA_MOTION_SRC} ${REL_JAZ_SPA_MOTION_H}
		${REL_JAZ_SPA_MATH_SRC} ${REL_JAZ_SPA_MATH_H}
		${REL_JAZ_SPA_IMG_PROC_SRC} ${REL_JAZ_SPA_IMG_PROC_H}
		${REL_JAZ_SPA_TOMO_SRC} ${REL_JAZ_SPA_TOMO_H}
		${REL_JAZ_SPA_IO_SRC} ${REL_JAZ_SPA_IO_H}
		${REL_JAZ_SPA_PROGRAMS_SRC} ${REL_JAZ_SPA_PROGRAMS_H}
		${REL_D3x3_SRC} ${REL_D3x3_H}
		${REL_LBFGS_SRC} ${REL_LBFGS_H} 
		${REL_SPHERICAL_HARMONICS_SRC} ${REL_SPHERICAL_HARMONICS_H})

	if(GUI)
		add_library(relion_gui_lib STATIC ${REL_GUI_SRC} ${REL_SRC_H} ${REL_HP})
	endif(GUI)
endif()

if(NOT MKLFFT)
	target_link_libraries(relion_lib ${FFTW_LIBRARIES})
	if(BUILD_OWN_FFTW)
		add_dependencies(relion_lib own_fftw_lib)
	endif()
	if(BUILD_OWN_FFTWF)
		add_dependencies(relion_lib own_fftwf_lib)
	endif()
endif(NOT MKLFFT)

if(GUI)
	include_directories("${FLTK_INCLUDE_DIR}")
	target_link_libraries(relion_gui_lib relion_lib ${FLTK_LIBRARIES})
	if(BUILD_OWN_FLTK)
		add_dependencies(relion_gui_lib OWN_FLTK)
	endif()
endif(GUI)

if (CUDA_FOUND)
	file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/acc/cuda/*.cu" "${CMAKE_SOURCE_DIR}/src/acc/cuda/cuda_kernels/*.cu" )
	cuda_add_library(relion_gpu_util ${REL_CUDA_SRC})

	if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "profiling")
		find_library(NV_TOOLS_LIBRARIES   NAMES nvToolsExt  PATHS ${CUDA_TOOLKIT_ROOT_DIR}/lib ${CUDA_TOOLKIT_ROOT_DIR}/lib64)
		list(APPEND EXTRA_LIBS "${NV_TOOLS_LIBRARIES}")
		target_link_libraries(relion_gpu_util ${NV_TOOLS_LIBRARIES})
		message(STATUS "Adding extra library for NVIDIA profiling: ${NV_TOOLS_LIBRARIES}")
	endif()

	# Presently we have a number of (bad) circular dependencies between the gpu util
	# and relion libraries, which cause errors at least on OS X with clang. Tell the
	# compiler to ignore them.
	if(APPLE)
		set(new_link_flags "-undefined suppress -flat_namespace")
		get_target_property(existing_link_flags relion_gpu_util LINK_FLAGS)
		if(existing_link_flags)
			set(new_link_flags "${existing_link_flags} ${new_link_flags}")
		endif()
		set_target_properties(relion_gpu_util PROPERTIES LINK_FLAGS "${new_link_flags}")
	endif()

	list(APPEND EXTRA_LIBS "${CUDA_CUFFT_LIBRARIES}")
	if(BUILD_SHARED_LIBS)
		install (TARGETS relion_gpu_util LIBRARY DESTINATION lib)
	else()
		target_link_libraries(relion_gpu_util relion_lib)
		target_link_libraries(relion_gpu_util ${CUDA_CUFFT_LIBRARIES})
	endif()

	target_link_libraries(relion_lib relion_gpu_util ${CUDA_CUFFT_LIBRARIES})
	target_link_libraries(relion_lib relion_gpu_util ${CUDA_CUFFT_LIBRARIES} ${CUDA_curand_LIBRARY})



	file(GLOB REL_JAZ_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/jaz/cuda/*.cu" "${CMAKE_SOURCE_DIR}/src/jaz/cuda/kernels/*.cu" )
	cuda_add_library(relion_jaz_gpu_util ${REL_JAZ_CUDA_SRC})

	#list(APPEND EXTRA_LIBS "${CUDA_CUFFT_LIBRARIES}")
	#if(BUILD_SHARED_LIBS)
	#	install (TARGETS relion_jaz_gpu_util LIBRARY DESTINATION lib)
	#else()
	#	target_link_libraries(relion_jaz_gpu_util relion_lib)
	#	target_link_libraries(relion_jaz_gpu_util ${CUDA_CUFFT_LIBRARIES})
	#endif()

	target_link_libraries(relion_lib relion_jaz_gpu_util ${CUDA_CUFFT_LIBRARIES})
	target_link_libraries(relion_lib relion_jaz_gpu_util ${CUDA_CUFFT_LIBRARIES} ${CUDA_curand_LIBRARY})


endif(CUDA_FOUND)

if(TIFF_FOUND)
	#message("TIFF FOUND")
	include_directories(${TIFF_INCLUDE_DIRS})
	target_link_libraries(relion_lib ${TIFF_LIBRARIES})
else()
	#message("TIFF NOT FOUND")
endif()

find_package(ZLIB)

find_package(PNG)
if(PNG_FOUND)
	add_definitions(-DHAVE_PNG)
	message("PNG FOUND")
	include_directories(${PNG_INCLUDE_DIRS})
	target_link_libraries(relion_lib ${PNG_LIBRARY})
else()
	message("PNG NOT FOUND")
endif()

find_package(JPEG)
if(JPEG_FOUND)
	add_definitions(-DHAVE_JPEG)
	message("JPEG FOUND")
	include_directories(${JPEG_INCLUDE_DIRS})
	target_link_libraries(relion_lib ${JPEG_LIBRARY})
else()
	message("JPEG NOT FOUND")
endif()

if(BUILD_OWN_TBB)
	add_dependencies(relion_lib OWN_TBB)
endif()

foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension

	if(${_target} STREQUAL "maingui") # We no longer want "relion_maingui"
		set(_target "relion")
		add_executable(${_target} maingui.cpp)
	else()
		add_executable(${_target} ${_target}.cpp )
		set_target_properties(${_target} PROPERTIES PREFIX "relion_")
	endif()

	set(LIB relion_lib)
	
 	add_dependencies(${_target} relion_lib)

	if(NOT MKLFFT)
		target_link_libraries(${_target} ${LIB} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	else()	
		target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	endif(NOT MKLFFT)

	if(CUDA_FOUND)
		target_link_libraries(${_target} relion_gpu_util)
	endif(CUDA_FOUND)
	
	if (ALTCPU)
		target_link_libraries(${_target} ${TBB_LIBRARIES})
	endif(ALTCPU)

	set_target_properties(${_target} PROPERTIES
	                      XX_STANDARD 11
	                      XX_STANDARD_REQUIRED ON
	                      XX_EXTENSIONS OFF
	)
	
	list(FIND GUI_TARGETS ${_target} IS_GUI_TARGET)
	if((NOT ${IS_GUI_TARGET} LESS 0) OR (${_target} STREQUAL "relion"))
		add_dependencies(${_target} relion_gui_lib)
		target_link_libraries(${_target} relion_gui_lib ${FLTK_LIBRARIES} ${X11})
	endif()

	if(TIFF_FOUND)
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()
    
	install (TARGETS ${_target} RUNTIME DESTINATION bin)

endforeach()

# Python script for running pytorch is required for class_ranker along with torch model
if(FETCH_TORCH_MODELS)
	add_dependencies(class_ranker copy_scripts class_ranker_model_file)
endif()

foreach (_target ${RELION_JAZ_SPA_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)	
	add_executable(${_target} ${CMAKE_SOURCE_DIR}/src/jaz/single_particle/apps/${_target}.cpp)
	set(LIB relion_lib)	
 	add_dependencies(${_target} relion_lib)	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_")
	target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	set_target_properties(${_target} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF	)
	
	if(TIFF_FOUND)
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()

	if (ALTCPU)
		target_link_libraries(${_target} ${TBB_LIBRARIES})
	endif(ALTCPU)
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

# former Dynamite programs:
foreach (_target ${RELION_JAZ_TOMO_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_tomo_${_target}" NAME_WE)	
	add_executable(${_target} ${CMAKE_SOURCE_DIR}/src/jaz/tomography/apps/${_target}.cpp)
	set(LIB relion_lib)	
 	add_dependencies(${_target} relion_lib)	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_tomo_")
	target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	set_target_properties(${_target} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF	)
	
	if(TIFF_FOUND)
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

foreach (_target ${RELION_JAZ_EXPERIMENTAL_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_exp_${_target}" NAME_WE)	
	add_executable(${_target} ${CMAKE_SOURCE_DIR}/src/jaz/scripts/${_target}.cpp)		
	set(LIB relion_lib)	
 	add_dependencies(${_target} relion_lib)	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_exp_")
	target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	set_target_properties(${_target} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON CXX_EXTENSIONS OFF	)
	
	if(TIFF_FOUND)
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()

	if (ALTCPU)
		target_link_libraries(${_target} ${TBB_LIBRARIES})
	endif(ALTCPU)
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	#message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	target_link_libraries(relion_lib ${OpenMP_omp_LIBRARY})
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# Set this flag to activate bounds checking in stl-vectors (incl. strings)
# It is useful to do this periodically, as it catches 
# difficult-to-see and rare-to-manifest bugs
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_DEBUG")
