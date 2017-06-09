include_directories("${CMAKE_SOURCE_DIR}")

include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/acc/*.cpp" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.cpp")
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/acc/*.h" "${CMAKE_SOURCE_DIR}/src/acc/cpu/*.h")
file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")

file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")

#--Remove apps using X11 if no GUI--
if(NOT GUI)
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/manualpicker.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_entries.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_jobwindow.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_mainwindow.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/displayer.cpp")
    list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/display.cpp")
    list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/maingui.cpp")
    list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/manualpick.cpp")
    list(REMOVE_ITEM RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/pipeliner.cpp")
endif(NOT GUI)


# relion_lib is STATIC or SHARED type based on BUILD_SHARED_LIBS=ON/OFF
if(BUILD_SHARED_LIBS)
    add_library(relion_lib SHARED ${REL_SRC} ${REL_SRC_H} ${REL_HP})
    install(TARGETS relion_lib LIBRARY DESTINATION lib)
else()
    add_library(relion_lib STATIC ${REL_SRC} ${REL_SRC_H} ${REL_HP})
endif()

target_link_libraries(relion_lib ${FFTW_LIBRARIES})
if(BUILD_OWN_FFTW)
    add_dependencies(relion_lib OWN_FFTW)
endif()

if(GUI)
    target_link_libraries(relion_lib ${FLTK_LIBRARIES})
    if(BUILD_OWN_FLTK)
        add_dependencies(relion_lib OWN_FLTK)
    endif()
endif(GUI)



if (CUDA_FOUND)

    file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/acc/cuda/*.cu" "${CMAKE_SOURCE_DIR}/src/acc/cuda/cuda_kernels/*.cu")
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
	
	if(GUI)
        include_directories("${FLTK_INCLUDE_DIR}")
		target_link_libraries(${_target} ${FLTK_LIBRARIES} ${X11})
	endif(GUI)
	
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


