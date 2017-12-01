include_directories("${CMAKE_SOURCE_DIR}")

include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cpp")
file(GLOB REL_GUI_SRC "${CMAKE_SOURCE_DIR}/src/manualpicker.cpp" "${CMAKE_SOURCE_DIR}/src/gui_*.cpp" "${CMAKE_SOURCE_DIR}/src/displayer.cpp")
# Remove GUI files from relion_lib
foreach(GUI_SRC_FILE ${REL_GUI_SRC})
	list(REMOVE_ITEM REL_SRC "${GUI_SRC_FILE}")
endforeach()
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.h")
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
    add_library(relion_lib SHARED ${REL_SRC} ${REL_SRC_H} ${REL_HP})
    install(TARGETS relion_lib LIBRARY DESTINATION lib)
    if(GUI)
        add_library(relion_gui_lib SHARED ${REL_GUI_SRC} ${REL_SRC_H} ${REL_HP})
        install(TARGETS relion_gui_lib LIBRARY DESTINATION lib)
    endif(GUI)
else()
    add_library(relion_lib STATIC ${REL_SRC} ${REL_SRC_H} ${REL_HP})
    if(GUI)
        add_library(relion_gui_lib STATIC ${REL_GUI_SRC} ${REL_SRC_H} ${REL_HP})
    endif(GUI)
endif()

target_link_libraries(relion_lib ${FFTW_LIBRARIES})
if(BUILD_OWN_FFTW)
    add_dependencies(relion_lib OWN_FFTW)
endif()

if(GUI)
    include_directories("${FLTK_INCLUDE_DIR}")
    target_link_libraries(relion_gui_lib relion_lib ${FLTK_LIBRARIES})
    if(BUILD_OWN_FLTK)
        add_dependencies(relion_gui_lib OWN_FLTK)
    endif()
endif(GUI)



if (CUDA_FOUND)

    file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cu" "${CMAKE_SOURCE_DIR}/src/gpu_utils/cuda_kernels/*.cu")
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
    
endif(CUDA_FOUND)

foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
	
    add_executable(${_target} ${_target}.cpp )
    set(LIB relion_lib)
	
 	add_dependencies(${_target} relion_lib)
	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_")
	
	target_link_libraries(${_target} ${LIB} ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	
	if(CUDA_FOUND)
	    target_link_libraries(${_target} relion_gpu_util)
	endif(CUDA_FOUND)
	
#	if(${_target} IN_LIST GUI_TARGETS)
#        add_dependencies(${_target} relion_gui_lib)
#		target_link_libraries(${_target} relion_gui_lib ${FLTK_LIBRARIES} ${X11})
#	endif()
	
	list(FIND GUI_TARGETS ${_target} IS_GUI_TARGET)
	if(NOT ${IS_GUI_TARGET} LESS 0)
        add_dependencies(${_target} relion_gui_lib)
        target_link_libraries(${_target} relion_gui_lib ${FLTK_LIBRARIES} ${X11})
    endif()
    
	message(STATUS "added ${_target}...")
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

if(GUI)
    ADD_CUSTOM_COMMAND(TARGET maingui
              POST_BUILD
              COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/bin/relion_maingui ${PROJECT_BINARY_DIR}/bin/relion
	      COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/bin/relion_qsub.csh ${PROJECT_BINARY_DIR}/bin/qsub.csh
    )
endif(GUI)



