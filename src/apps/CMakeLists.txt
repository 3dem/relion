include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

#set(CUDA_NVCC_FLAGS "-g -G")

find_path (X11_INCLUDES Xdbe.h)

#FLTK_WRAP_UI(Ui ui.fl)
#ADD_LIBRARY(Ui callbacks.cxx ${Ui_FLTK_UI_SRCS})

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cu")
file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cpp")
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.h")
file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")
file(GLOB REL_LIB "${CMAKE_SOURCE_DIR}/lib/*.so" )

file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")

add_library(relion_lib OBJECT ${REL_SRC} ${REL_SRC_H} ${REL_HP})
foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
	if(CUDA_FOUND)
	    cuda_add_executable(${_target} ${_target}.cpp ${REL_CUDA_SRC} $<TARGET_OBJECTS:relion_lib>)
	else(CUDA_FOUND)
	    add_executable(${_target} ${_target}.cpp $<TARGET_OBJECTS:relion_lib>)
	endif(CUDA_FOUND)
	#add_dependencies(${_target} FFTW Ui)
	#target_link_libraries(${_target} Ui ${FLTK_LIBRARIES} ${FFTW_LIBRARIES} ${X11} ${relion_lib} )
	target_link_libraries(${_target} ${FLTK_LIBRARIES} ${FFTW_LIBRARIES} ${X11} ${relion_lib} )
	message(STATUS "added ${_target}...")
endforeach()
