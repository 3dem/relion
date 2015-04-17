include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${FFTW_INCLUDES}")
#include_directories(${CMAKE_BINARY_DIR}/include)

find_path (X11_INCLUDES Xdbe.h)

#FLTK_WRAP_UI(Ui ui.fl)
#ADD_LIBRARY(Ui callbacks.cxx ${Ui_FLTK_UI_SRCS})

message(STATUS "CMAKE_BINARY_DIR:" ${CMAKE_BINARY_DIR})

file(GLOB REL_SRC "${CMAKE_SOURCE_DIR}/src/*.cpp" "${CMAKE_SOURCE_DIR}/src/*.c" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cpp")
file(GLOB REL_SRC_H "${CMAKE_SOURCE_DIR}/src/*.h" "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.h")
file(GLOB REL_HP "${CMAKE_SOURCE_DIR}/src/Healpix_2.15a/*.cc")

#--Apps using X11--
if(BUILD_GUI_TARGETS)
    file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/refine.cpp" "${CMAKE_SOURCE_DIR}/src/apps/display.cpp")
    #file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/*.cpp")
else(BUILD_GUI_TARGETS)
    file(GLOB RELION_TARGETS "${CMAKE_SOURCE_DIR}/src/apps/refine.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/manualpicker.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_entries.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_jobwindow.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/gui_mainwindow.cpp")
    list(REMOVE_ITEM REL_SRC "${CMAKE_SOURCE_DIR}/src/displayer.cpp")
endif(BUILD_GUI_TARGETS)

if (CUDA_FOUND)
    file(GLOB REL_CUDA_SRC "${CMAKE_SOURCE_DIR}/src/gpu_utils/*.cu")
    cuda_add_library(gpu_util STATIC ${REL_CUDA_SRC})
endif(CUDA_FOUND)

add_library(relion_lib OBJECT ${REL_SRC} ${REL_SRC_H} ${REL_HP})

foreach (_target ${RELION_TARGETS})
	GET_FILENAME_COMPONENT(_target "relion_${_target}" NAME_WE)  #specify target name WE=WithoutExtension
    add_executable(${_target} ${_target}.cpp $<TARGET_OBJECTS:relion_lib>)
	
	if(BUILD_GUI_TARGETS)
	    target_link_libraries(${_target} ${FLTK_LIBRARIES} gpu_util ${FFTW_LIBRARIES} ${X11} ${relion_lib} )
	else(BUILD_GUI_TARGETS)
	    target_link_libraries(${_target} gpu_util ${FFTW_LIBRARIES} ${relion_lib} )
	endif(BUILD_GUI_TARGETS)
	
	message(STATUS "added ${_target}...")
endforeach()
