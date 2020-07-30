#################################################################################
# Description:  Cmake helper to download and configure LibTorch at cmake config
#               time. This is done WITHOUT the cmake util ExternalProject,
#               because LibTorch is required during config time. It is instead
#               implemented using linux command execution and is thus not
#               cross-platform compatible. Download integrity is managed with
#               MD5 checksum and 5 attempts will be made to download and unpack.
# Author:       Dari Kimanius (github.com/dkimanius)
# Date created: 2020-06-07
#################################################################################

set(TORCH_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/libtorch")
find_package(Torch PATHS ${TORCH_EXTERNAL_PATH})
option(TORCH_EXTERNAL_VERBOSE "Verbose mode during Torch fetch and configure." OFF)
if(NOT TORCH_FOUND)
	set(TORCH_FILE_NAME libtorch-commit-4c7fb8c2-cxx11-abi-glibc2.17-shared.zip)

	set(TORCH_URL "ftp://ftp.mrc-lmb.cam.ac.uk/pub/dari/${TORCH_FILE_NAME}")
	set(TORCH_HASH 4aa641cd867d141e438596a1ba31e13d)

        set(TORCH_FILE_PATH ${TORCH_EXTERNAL_PATH}/${TORCH_FILE_NAME})
	
	message(STATUS "Downloading Torch...")
	file(MAKE_DIRECTORY ${TORCH_EXTERNAL_PATH})
	foreach(ATTEMPT RANGE 2 6)
		if (BUILD_TORCH_VERBOSE)
			execute_process(
					COMMAND wget -O ${TORCH_FILE_PATH} ${TORCH_URL}
					RESULT_VARIABLE WGET_FAIL)
		else()
			execute_process(
					COMMAND wget -O ${TORCH_FILE_PATH} ${TORCH_URL}
					RESULT_VARIABLE WGET_FAIL
					OUTPUT_QUIET ERROR_QUIET)
		endif()
		if(NOT WGET_FAIL)
			file(MD5 ${TORCH_FILE_PATH} CHECKSUM)
			if (${CHECKSUM} STREQUAL ${TORCH_HASH})
				if (BUILD_TORCH_VERBOSE)
					execute_process(
							COMMAND unzip libtorch/${TORCH_FILE_NAME}
							WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external
							RESULT_VARIABLE UNZIP_FAIL)
				else()
					execute_process(
							COMMAND unzip libtorch/${TORCH_FILE_NAME}
							WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/external
							RESULT_VARIABLE UNZIP_FAIL
							OUTPUT_QUIET ERROR_QUIET)
				endif()
				if(NOT UNZIP_FAIL)
					message(STATUS "Torch was successfully downloaded and configured")
					break()
				else()
					message(FATAL_ERROR "Unpacking failed with error: ${UNZIP_FAIL}")
				endif()
			else()
				message(STATUS "Torch checksum mismatch")
			endif()
		else()
			message(STATUS "Download failed with error: ${WGET_FAIL}")
		endif()
		if (ATTEMPT LESS 6)
			message(STATUS "Retrying (attempt ${ATTEMPT})...")
		else()
			message(FATAL_ERROR "Could not donwload and configure TORCH")
		endif()
	endforeach()
	find_package(Torch REQUIRED PATHS ${TORCH_EXTERNAL_PATH})
else(NOT TORCH_FOUND)
	message(STATUS "Found previously built non-system Torch libraries that will be used.")
endif()

include_directories("${TORCH_INCLUDE_DIRS}")
link_directories(${TORCH_LIBRARY_DIRS})

install(DIRECTORY ${TORCH_LIBRARY_DIRS}/ DESTINATION lib
	USE_SOURCE_PERMISSIONS FILES_MATCHING PATTERN "*.so*")

#message(STATUS "TORCH_INCLUDE_DIRS: ${TORCH_INCLUDE_DIRS}")
#message(STATUS "TORCH_LIBRARY_DIRS: ${TORCH_LIBRARY_DIRS}")
