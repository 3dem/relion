
set(TBB_PREFIX tbb2020U3)
set(TBB_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/tbb")

## ------------------------------------------------------------- PREVIOUS EXT LIBS? --

find_library(TBB_TEST_LIB NAMES tbb PATHS "${TBB_EXTERNAL_PATH}/${TBB_PREFIX}/build/${TBB_PREFIX}_release" NO_DEFAULT_PATH)
find_path(TBB_TEST_INCLUDES NAMES tbb/tbb.h PATHS "${TBB_EXTERNAL_PATH}/${TBB_PREFIX}/include" NO_DEFAULT_PATH)

if(TBB_TEST_LIB AND TBB_TEST_INCLUDES)
    set(TBB_FOUND TRUE)
	message(STATUS "Found previously built non-system TBB libraries that will be used.")
else()
	message(STATUS "TBB_TEST_LIB: ${TBB_TEST_LIB}")
	message(STATUS "TBB_TEST_INCLUDES: ${TBB_TEST_INCLUDES}")
    set(TBB_FOUND FALSE)
	message(STATUS "--------------------------------------------------------")
	message(STATUS "-------- NO EXISTING TBB LIBRARIES WHERE FOUND. -------")
	message(STATUS "-------------- TBB WILL BE DOWNLOADED AND -------------")
	message(STATUS "--------------- BUILT DURING COMPILE-TIME. -------------")
	message(STATUS "--------------------------------------------------------")
	message(STATUS "---- A WORKING INTERNET CONNECTION WILL BE REQUIRED. ---")
	message(STATUS "--------------------------------------------------------")
endif()

## ----------------------------------------------------------------- NEW EXT LIBS? --  
 
if(NOT TBB_FOUND)
    message(STATUS "no previous tbb found, the following paths are set for libs/headers TO BE built")
    
    include(ExternalProject)

    set(TBB_URL https://github.com/oneapi-src/oneTBB/archive/v2020.3.tar.gz)   # TBB 2020 U3
    set(TBB_URL_MD5 ea8fa4332f4bad10a75a361cba025380)
    set(TBB_TAR_NAME tbb-2020_U3.tar.gz)
   
    ExternalProject_Add(OWN_TBB
		URL ${TBB_URL}
		URL_MD5 ${TBB_URL_MD5}
		DOWNLOAD_DIR ${TBB_EXTERNAL_PATH}
		DOWNLOAD_NAME ${TBB_TAR_NAME}
		SOURCE_DIR "${TBB_EXTERNAL_PATH}/${TBB_PREFIX}"
		CONFIGURE_COMMAND ""
		BUILD_COMMAND make tbb_build_prefix=${TBB_PREFIX}
		BUILD_IN_SOURCE 1
		INSTALL_COMMAND ""
		LOG_DOWNLOAD 1
		LOG_BUILD 1
	)
	
    set(BUILD_OWN_TBB TRUE)
else(NOT TBB_FOUND)
    set(BUILD_OWN_TBB FALSE)
endif(NOT TBB_FOUND)

set(TBB_INCLUDE_DIRS "${TBB_EXTERNAL_PATH}/${TBB_PREFIX}/include" )
include_directories("${TBB_INCLUDE_DIRS}")
	
# in release mode
set(TBB_LIBRARY_DIRS ${TBB_EXTERNAL_PATH}/${TBB_PREFIX}/build/${TBB_PREFIX}_release)
link_directories(${TBB_LIBRARY_DIRS})
set(TBB_LIBRARIES tbb tbbmalloc)

# in debug mode
#set(TBB_LIBRARY_DIRS ${TBB_EXTERNAL_PATH}/build/${TBB_PREFIX}_debug)
#link_directories(${TBB_LIBRARY_DIR})
#set(TBB_LIBRARIES tbb_debug tbbmalloc_debug)
	
install(DIRECTORY ${TBB_LIBRARY_DIRS}/ DESTINATION lib
	USE_SOURCE_PERMISSIONS FILES_MATCHING PATTERN "*.so*")

message(STATUS "TBB_INCLUDE_DIRS: ${TBB_INCLUDE_DIRS}")
message(STATUS "TBB_LIBRARY_DIRS: ${TBB_LIBRARY_DIRS}")
