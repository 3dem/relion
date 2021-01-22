set(FFTW_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fftw")

if (NOT DEFINED TARGET_X86)
    try_compile(TARGET_X86 ${CMAKE_BINARY_DIR}
                "${CMAKE_SOURCE_DIR}/cmake/TestX86.c")
endif()

find_path(   OWN_FFTW_INCLUDES NAMES fftw3.h PATHS ${FFTW_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
find_library(OWN_FFTW_SINGLE   NAMES fftw3f  PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)
find_library(OWN_FFTW_DOUBLE   NAMES fftw3   PATHS ${FFTW_EXTERNAL_PATH}/lib     NO_DEFAULT_PATH)

if(OWN_FFTW_INCLUDES AND (OWN_FFTW_SINGLE OR NOT FFTW_SINGLE_REQUIRED) AND (OWN_FFTW_DOUBLE OR NOT FFTW_DOUBLE_REQUIRED))

	if (OWN_FFTW_SINGLE AND FFTW_SINGLE_REQUIRED)
		message(STATUS "Found previously built non-system single precision FFTW libraries that will be used.")
		#message(STATUS "OWN_FFTW_SINGLE:    ${OWN_FFTW_SINGLE}")
	endif()

	if (OWN_FFTW_DOUBLE AND FFTW_DOUBLE_REQUIRED)
		message(STATUS "Found previously built non-system double precision FFTW libraries that will be used.")
		#message(STATUS "OWN_FFTW_DOUBLE:    ${OWN_FFTW_DOUBLE}")
	endif()
	
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
		if (AMDFFTW)
			set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-sse2 --enable-avx --enable-avx2 --enable-amd-opt)
		else()
			set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-avx)
		endif()
	endif()
	
	set(OWN_FFTW_SINGLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_DOUBLE ${FFTW_EXTERNAL_PATH}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
	set(OWN_FFTW_INCLUDES "${FFTW_EXTERNAL_PATH}/include" )
	
	set(FFTW_PATH ${FFTW_PATH} ${FFTW_EXTERNAL_PATH})

	set(FFTW_EXTERNAL_LIBS_TAR_DIRECTORY  ${FFTW_EXTERNAL_PATH})
	set(FFTW_EXTERNAL_LIBS_EXTRACT_TARGET ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY})

	set(FFTW_FFTW3_TAR_FILE http://fftw.org/fftw-3.3.8.tar.gz)
	set(FFTW_MD5 8aac833c943d8e90d51b697b27d4384d)

	if (AMDFFTW)
		set(FFTW_FFTW3_TAR_FILE https://github.com/amd/amd-fftw/archive/2.2.zip)
		set(FFTW_MD5 2e9c59ad80ec5bd75ce04c7970c9f47a)
	endif()

	set(FFTW_FFTW3_LIB_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3)
	set(FFTW_FFTW3_BUILD_DIR ${FFTW_EXTERNAL_LIBS_EXTRACT_TARGET}/fftw3-build)

	include(ExternalProject)
	
	if (FFTW_SINGLE_REQUIRED)
		set(BUILD_OWN_FFTW TRUE)
	endif()
	
	if (FFTW_DOUBLE_REQUIRED)
		set(BUILD_OWN_FFTWF TRUE)
	endif()

	# Rather messy logic:
	# We build double prec here but if double prec is not required, build single prec	
	if (NOT FFTW_DOUBLE_REQUIRED)
		set(ext_conf_flags_fft ${ext_conf_flags_fft} --enable-float --enable-sse)
	endif()

	externalproject_add(own_fftw_lib
	URL ${FFTW_FFTW3_TAR_FILE}
	URL_MD5 ${FFTW_MD5}
	DOWNLOAD_DIR ${FFTW_EXTERNAL_LIBS_TAR_DIRECTORY}
	SOURCE_DIR ${FFTW_FFTW3_LIB_DIR}
	CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fft}
	INSTALL_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BINARY_DIR ${FFTW_EXTERNAL_PATH}/fftw3
	BUILD_COMMAND ${MAKE}
	#	LOG_CONFIGURE
	#	LOG_BUILD
	LOG_INSTALL)

	add_custom_command(
		COMMAND ${CMAKE_COMMAND} -E echo "Registering own FFTW byproducts"
		OUTPUT "${FFTW_EXTERNAL_PATH}/lib/libfftw3.so"
		DEPENDS own_fftw_lib)
	add_custom_target(own_fftw_lib_byproducts
		DEPENDS "${FFTW_EXTERNAL_PATH}/lib/libfftw3.so")

	# When both double and single prec are required, build single later.	
	if (FFTW_DOUBLE_REQUIRED AND FFTW_SINGLE_REQUIRED)
		externalproject_add(own_fftwf_lib
		URL ${FFTW_FFTW3_TAR_FILE}
		URL_MD5 ${FFTW_MD5}
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

		add_custom_command(
			COMMAND ${CMAKE_COMMAND} -E echo "Registering own FFTWf byproducts"
			OUTPUT "${FFTW_EXTERNAL_PATH}/lib/libfftw3f.so"
			DEPENDS own_fftwf_lib)
		add_custom_target(own_fftwf_lib_byproducts
			DEPENDS "${FFTW_EXTERNAL_PATH}/lib/libfftw3f.so")
	endif()

endif()

if (FFTW_SINGLE_REQUIRED)
	set(FFTW_LIBRARIES ${OWN_FFTW_SINGLE} ${FFTW_LIBRARIES})
endif()

if (FFTW_DOUBLE_REQUIRED)
	set(FFTW_LIBRARIES ${OWN_FFTW_DOUBLE} ${FFTW_LIBRARIES})
endif()

if (FFTW_INCLUDES)
	set(FFTW_INCLUDES ${OWN_FFTW_INCLUDES} ${FFTW_INCLUDES})
else()
	set(FFTW_INCLUDES ${OWN_FFTW_INCLUDES})
endif()

#message(STATUS "FFTW_INCLUDES:     ${FFTW_INCLUDES}")
#message(STATUS "FFTW_LIBRARIES:    ${FFTW_LIBRARIES}")
