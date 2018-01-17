set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")

if (NOT DEFINED TARGET_X86)
    try_compile(TARGET_X86 ${CMAKE_BINARY_DIR}
                "${CMAKE_SOURCE_DIR}/cmake/TestX86.c")
endif()

find_path(   OWN_FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
find_library(OWN_FFTW_SINGLE   NAMES fftw3   PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)
find_library(OWN_FFTW_DOUBLE   NAMES fftw3f  PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)

if(OWN_FFTW_INCLUDES AND (OWN_FFTW_SINGLE OR NOT FFTW_SINGLE_REQUIRED) AND (OWN_FFTW_DOUBLE OR NOT FFTW_DOUBLE_REQUIRED))

	message(STATUS "Found previously built non-system FFTW libraries that will be used.")

	set(FFTW_FOUND TRUE)
	set(BUILD_OWN_FFTW FALSE)
	set(BUILD_OWN_FFTWF FALSE)
	
else()
	
	message(STATUS "--------------------------------------------------------")
	message(STATUS "------ REQUIRED FFTW LIBRARIES WHERE NOT FOUND. --------")
	message(STATUS "-------------- FFTW WILL BE DOWNLOADED AND -------------")
	message(STATUS "--------------- BUILT DURING COMPILE-TIME. -------------")
	message(STATUS "--------------------------------------------------------")
	message(STATUS "---- A WORKING INTERNET CONNECTION WILL BE REQUIRED. ---")
	message(STATUS "--------------------------------------------------------")
	
	set(FFTW_FOUND FALSE)
	
	set(ext_conf_flags_fft --enable-shared --prefix=${FFTW_EXTERNAL_PATH})
	if(TARGET_X86)
		set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-avx)
	endif()
	
	set(OWN_FFTW_SINGLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_DOUBLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_INCLUDES "${FFTW_EXTERNAL_PATH}/include" )
	
	set(FFTW_PATH ${FFTW_PATH} ${FFTW_EXTERNAL_PATH})

	set(FFTW_EXTERNAL_LIBS_TAR_DIRECTORY  ${FFTW_EXTERNAL_PATH})
	set(FFTW_EXTERNAL_LIBS_EXTRACT_TARGET ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY})

	set(FFTW_FFTW3_TAR_FILE ftp://ftp.fftw.org/pub/fftw/fftw-3.3.4.tar.gz)
	#set(FFTW_FFTW3_TAR_FILE https://drive.google.com/uc?export=download&id=0B942d76zVnSeazZWcExRaXIyVDg) #backup location fftw-3.3.4
	#set(FFTW_TAR_NAME fftw-3.3.4.tar.gz)

	set(FFTW_FFTW3_LIB_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
	set(FFTW_FFTW3_BUILD_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)

	include(ExternalProject)
	
	if (FFTW_SINGLE_REQUIRED)
		set(BUILD_OWN_FFTW TRUE)
	endif()
	
	if (FFTW_DOUBLE_REQUIRED)
		set(BUILD_OWN_FFTWF TRUE)
	endif()
	
	if (NOT FFTW_DOUBLE_REQUIRED)
		set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-float --enable-sse)
	endif()
	
	externalproject_add(own_fftw_lib
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
	
	
	if (FFTW_DOUBLE_REQUIRED AND FFTW_SINGLE_REQUIRED)
		
		externalproject_add(own_fftwf_lib
		URL ${FFTW_FFTW3_TAR_FILE}
		URL_MD5 2edab8c06b24feeb3b82bbb3ebf3e7b3
		DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
		SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
		CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}  --enable-float --enable-sse
		INSTALL_DIR ${FFTW_EXTERNAL_PATH}/fftw3
		BINARY_DIR ${FFTW_EXTERNAL_PATH}/fftw3
		BUILD_COMMAND ${MAKE}
		#	LOG_CONFIGURE
		#	LOG_BUILD
		LOG_INSTALL)
	
		add_dependencies(own_fftwf_lib own_fftw_lib)
	endif()

endif()

if (FFTW_SINGLE_REQUIRED)
	set(FFTW_LIBRARIES ${OWN_FFTW_SINGLE} ${FFTW_LIBRARIES})
endif()

if (FFTW_DOUBLE_REQUIRED)
	set(FFTW_LIBRARIES ${OWN_FFTW_DOUBLE} ${FFTW_LIBRARIES})
endif()

set(FFTW_INCLUDES ${OWN_FFTW_INCLUDES} ${FFTW_INCLUDES})

#message(STATUS "FFTW_INCLUDES:     ${FFTW_INCLUDES}")
#message(STATUS "FFTW_LIBRARIES:    ${FFTW_LIBRARIES}")
