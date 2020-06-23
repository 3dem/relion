find_package(Catch2)

include(Catch)

add_executable(catchTest
	test_ctf.cpp
	)

add_dependencies(catchTest relion_lib)
target_link_libraries(catchTest Catch2::Catch2)
target_link_libraries(catchTest relion_lib)
target_link_libraries(catchTest ${FFTW_LIBRARIES})
target_link_libraries(catchTest ${TIFF_LIBRARIES})

include_directories(catchTest PRIVATE ${CMAKE_SOURCE_DIR})
include_directories(catchTest PRIVATE ${FFTW_INCLUDES})

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

catch_discover_tests(catchTest)
