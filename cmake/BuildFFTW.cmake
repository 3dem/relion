set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")

if (NOT DEFINED TARGET_X86)
    try_compile(TARGET_X86 ${CMAKE_BINARY_DIR}
                "${CMAKE_SOURCE_DIR}/cmake/TestX86.c")
endif()

find_path(   FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
find_library(FFTW_SINGLE   NAMES fftw3   PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)
find_library(FFTW_DOUBLE   NAMES fftw3f  PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)

if(FFTW_INCLUDES AND (FFTW_SINGLE OR NOT FFTW_SINGLE_REQUIRED) AND (FFTW_DOUBLE OR NOT FFTW_DOUBLE_REQUIRED))

	message(STATUS "Found previously built non-system FFTW libraries that will be used.")

	set(FFTW_FOUND TRUE)
	set(BUILD_OWN_FFTW FALSE)
	set(FFTW_LIBRARIES ${FFTW_SINGLE} ${FFTW_DOUBLE})
	
else()
	
	message(STATUS "--------------------------------------------------------")
	message(STATUS "-------- NO EXISTING FFTW LIBRARIES WHERE FOUND. -------")
	message(STATUS "-------------- FFTW WILL BE DOWNLOADED AND -------------")
	message(STATUS "--------------- BUILT DURING COMPILE-TIME. -------------")
	message(STATUS "--------------------------------------------------------")
	message(STATUS "---- A WORKING INTERNET CONNECTION WILL BE REQUIRED. ---")
	message(STATUS "--------------------------------------------------------")
	
	set(FFTW_FOUND FALSE)
	set(BUILD_OWN_FFTW TRUE)
	
	set(ext_conf_flags_fft --enable-shared --enable-float --prefix=${FFTW_EXTERNAL_PATH})
	if(TARGET_X86)
		set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-sse --enable-avx)
	endif()
	
	set(FFTW_SINGLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(FFTW_DOUBLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(FFTW_LIBRARIES ${FFTW_SINGLE} ${FFTW_DOUBLE})
	
	set(FFTW_PATH      "${FFTW_EXTERNAL_PATH}" )
	set(FFTW_INCLUDES  "${FFTW_EXTERNAL_PATH}/include" )

	include(ExternalProject)
	set(FFTW_EXTERNAL_LIBS_TAR_DIRECTORY  ${FFTW_EXTERNAL_PATH})
	set(FFTW_EXTERNAL_LIBS_EXTRACT_TARGET ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY})
	message( STATUS "no previous fftw found, the following paths are set for libs/headers TO BE built")

	set(FFTW_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
	#set(FFTW_FFTW3_TAR_FILE https://drive.google.com/uc?export=download&id=0B942d76zVnSeazZWcExRaXIyVDg) #backup location fftw-3.3.4
	#set(FFTW_TAR_NAME fftw-3.3.4.tar.gz)

	set(FFTW_FFTW3_LIB_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
	set(FFTW_FFTW3_BUILD_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)

	#set(CMAKE_INSTALL_PREFIX  ${FFTW_EXTERNAL_PATH})
	externalproject_add(OWN_FFTW
	URL ${FFTW_FFTW3_TAR_FILE}
	URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
	DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
	SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}
	INSTALL_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BINARY_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BUILD_COMMAND ${MAKE}
	#	LOG_CONFIGURE
	#	LOG_BUILD
	LOG_INSTALL)

endif()

message(STATUS "FFTW_INCLUDES:     ${FFTW_INCLUDES}")
message(STATUS "FFTW_LIBRARIES:    ${FFTW_LIBRARIES}")
