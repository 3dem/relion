include_directories("${CMAKE_SOURCE_DIR}")

message("Running apps/CMakeLists.txt...")

if(NOT MKLFFT)
    include_directories("${FFTW_INCLUDES}")
endif(NOT MKLFFT)
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_GUI_SRC "${CMAKE_SOURCE_DIR}/src/manualpicker.cpp" "${CMAKE_SOURCE_DIR}/src/gui_*.cpp" "${CMAKE_SOURCE_DIR}/src/displayer.cpp")

if (ALTCPU)
    file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.cpp" "${CMAKE_SOURCE_DIR}/src/acc/cpu/cpu_kernels/*.cpp" )
    file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/acc/*.h" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.h" "${CMAKE_SOURCE_DIR}/src/acc/cpu/cpu_kernels/*.h" ) 
else()
    file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/acc/*.cpp" )
    file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/acc/*.h" ) 
endif(ALTCPU)

# Remove GUI files from relion_lib
foreach(GUI_SRC_FILE ${REL_GUI_SRC})
	list(REMOVE_ITEM REL_SRC "${GUI_SRC_FILE}")
endforeach()

file(GLOB REL_JAZ_SRC "${CMAKE_SOURCE_DIR}/src/jaz/*.cpp")
file(GLOB REL_JAZ_SRC_H "${CMAKE_SOURCE_DIR}/src/jaz/*.h")
file(GLOB REL_D3x3_SRC "${CMAKE_SOURCE_DIR}/src/jaz/d3x3/*.c")
file(GLOB REL_D3x3_SRC_H "${CMAKE_SOURCE_DIR}/src/jaz/d3x3/*.h")

file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")
file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")

set(GUI_TARGETS display maingui manualpick pipeliner)

#--Remove apps using X11 if no GUI--
if(NOT GUI)
    foreach(TARGET ${GUI_TARGETS})
        list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/${TARGET}.cpp")
    endforeach()
endif(NOT GUI)

# relion_lib is STATIC or SHARED type based on BUILD_SHARED_LIBS=ON/OFF
# relion_lib only contains non-X11 parts
# relion_gui_lib is where the X11 code is placed
if(BUILD_SHARED_LIBS)
	add_library(relion_lib SHARED ${REL_SRC} ${REL_SRC_H} ${REL_HP} ${REL_JAZ_SRC} ${REL_JAZ_H} ${REL_D3x3_SRC} ${REL_D3x3_H})
	install(TARGETS relion_lib LIBRARY DESTINATION lib)
	if(GUI)
		add_library(relion_gui_lib SHARED ${REL_GUI_SRC} ${REL_SRC_H} ${REL_HP})
		install(TARGETS relion_gui_lib LIBRARY DESTINATION lib)
	endif(GUI)
else()
	add_library(relion_lib STATIC ${REL_SRC} ${REL_SRC_H} ${REL_HP} ${REL_JAZ_SRC} ${REL_JAZ_H} ${REL_D3x3_SRC} ${REL_D3x3_H})
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
    
endif(CUDA_FOUND)

if(TIFF_FOUND)
	message("TIFF FOUND")
    include_directories(${TIFF_INCLUDE_DIRS})
else()
	message("TIFF NOT FOUND")
endif()

foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
	
	add_executable(${_target} ${_target}.cpp )
	set(LIB relion_lib)
	
 	add_dependencies(${_target} relion_lib)
	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_")

	if(NOT MKLFFT)
	    target_link_libraries(${_target} ${LIB} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	else()	
	    target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	endif(NOT MKLFFT)
	
	if(TIFF_FOUND)
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()

	if(CUDA_FOUND)
	    target_link_libraries(${_target} relion_gpu_util)
	endif(CUDA_FOUND)
	
	if (ALTCPU)
		target_link_libraries(${_target} ${TBB_LIBRARIES})
	endif(ALTCPU)

	set_target_properties(${_target} PROPERTIES
		CXX_STANDARD 11
		CXX_STANDARD_REQUIRED ON
		CXX_EXTENSIONS OFF
	)
	
	list(FIND GUI_TARGETS ${_target} IS_GUI_TARGET)
	if(NOT ${IS_GUI_TARGET} LESS 0)
		add_dependencies(${_target} relion_gui_lib)
		target_link_libraries(${_target} relion_gui_lib ${FLTK_LIBRARIES} ${X11})
	endif()

	if(TIFF_FOUND)
		message("linking ${_target} to TIFF (${TIFF_LIBRARIES})")
		target_link_libraries(${_target} ${TIFF_LIBRARIES})
	endif()
    
	#message(STATUS "added ${_target}...")
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

if(GUI)
    ADD_CUSTOM_COMMAND(TARGET maingui
              POST_BUILD
              COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/bin/relion_maingui ${PROJECT_BINARY_DIR}/bin/relion
	      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/bin/relion_qsub.csh ${PROJECT_BINARY_DIR}/bin/qsub.csh
    )
endif(GUI)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
