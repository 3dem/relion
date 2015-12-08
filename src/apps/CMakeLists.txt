include_directories("${CMAKE_SOURCE_DIR}")

include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cpp")
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.h")
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
endif(NOT GUI)

if (CUDA_FOUND)
    file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cu" "${CMAKE_SOURCE_DIR}/src/gpu_utils/cuda_kernels/*.cu")
    cuda_add_library(relion_gpu_util ${REL_CUDA_SRC})
    
    if (${CMAKE_BUILD_TYPE_LOWER} STREQUAL "profiling")
        find_library(NV_TOOLS_LIBRARIES   NAMES nvToolsExt  PATHS $ENV{CUDA_HOME}/lib64)
        list(APPEND EXTRA_LIBS "${NV_TOOLS_LIBRARIES}")
        message(STATUS "Adding extra library for NVIDIA profiling: ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so")
    endif()
    
    find_library(CUFFT_LIBRARIES   NAMES cufft  PATHS $ENV{CUDA_HOME}/lib64)
    list(APPEND EXTRA_LIBS "${CUFFT_LIBRARIES}")
    if(BUILD_SHARED_LIBS)
        install (TARGETS relion_gpu_util LIBRARY DESTINATION lib)
    endif()
endif(CUDA_FOUND)

# relion_lib is STATIC or SHARED type based on BUILD_SHARED_LIBS=ON/OFF
if(BUILD_SHARED_LIBS)
    add_library(relion_lib SHARED ${REL_SRC} ${REL_SRC_H} ${REL_HP})
    if(GUI)
    	target_link_libraries(relion_lib ${FLTK_LIBRARIES})
    endif(GUI)
    install (TARGETS relion_lib LIBRARY DESTINATION lib)
else()
    add_library(relion_lib OBJECT ${REL_SRC} ${REL_SRC_H} ${REL_HP})
endif()


foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
	
    if(BUILD_SHARED_LIBS)
	    add_executable(${_target} ${_target}.cpp )
	    set(LIB relion_lib)
	else()
	    add_executable(${_target} ${_target}.cpp $<TARGET_OBJECTS:relion_lib>)
	    set(LIB ${relion_lib})
	endif()
	
	if(NEW_OWN_FFTW)
 	   add_dependencies(${_target} FFTW3)
 	endif()
	if(NEW_OWN_FLTK)
 	   add_dependencies(${_target} FLTK)
 	endif()
	
	set_target_properties(${_target} PROPERTIES PREFIX "relion_")
	if(GUI)
	       include_directories("${FLTK_INCLUDE_DIR}") 
           if (CUDA_FOUND)
  	       target_link_libraries(${_target} ${LIB} ${FLTK_LIBRARIES} relion_gpu_util ${FFTW_LIBRARIES} ${X11} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	    else(CUDA_FOUND)
  	       target_link_libraries(${_target} ${LIB} ${FLTK_LIBRARIES}          ${FFTW_LIBRARIES} ${X11} ${EXTRA_LIBS} ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})
	    endif(CUDA_FOUND)
	else(GUI)
        if (CUDA_FOUND)
	       target_link_libraries(${_target} ${LIB} relion_gpu_util ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES})
	    else(CUDA_FOUND)
	       target_link_libraries(${_target} ${LIB}          ${FFTW_LIBRARIES} ${EXTRA_LIBS} ${MPI_LIBRARIES})
	    endif(CUDA_FOUND)
	endif(GUI)
	message(STATUS "added ${_target}...")
	
	install (TARGETS ${_target} RUNTIME DESTINATION bin)
endforeach()

if(GUI)
    ADD_CUSTOM_COMMAND(TARGET maingui
              POST_BUILD
              COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/bin/relion_maingui ${PROJECT_BINARY_DIR}/bin/relion
    )
endif(GUI)


