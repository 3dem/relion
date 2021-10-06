#################################################################################
# Description:  Cmake helper to download and configure Torch models.
# Author:       Dari Kimanius (github.com/dkimanius)
# Date created: 2021-05-04
#################################################################################

include(ExternalProject)
set(CLASS_RANKER_DEFAULT_MODEL_FILE_NAME "relion_class_ranker_default_model.pt")

set(MODELS_DIR  ${CMAKE_SOURCE_DIR}/external/torch_models)
if(NOT EXISTS ${MODELS_DIR})
	file(MAKE_DIRECTORY ${MODELS_DIR})
endif()

set(CLASS_RANKER_MODEL_FILE_NAME_TAR class_ranker_0.1.3_torch_1.0.1.pt.tar.gz)
set(CLASS_RANKER_MODEL_FILE_NAME class_ranker_0.1.3_torch_1.0.1.pt)
set(CLASS_RANKER_MODEL_URL "ftp://ftp.mrc-lmb.cam.ac.uk/pub/dari/${CLASS_RANKER_MODEL_FILE_NAME_TAR}")
set(CLASS_RANKER_MODEL_MD5 b39f0cbc31f510e3a03e89f1e11110fe)

set(CLASS_RANKER_MODEL_FOUND 0)

message(STATUS "Checking class ranker model file...")
if(EXISTS ${MODELS_DIR}/${CLASS_RANKER_MODEL_FILE_NAME_TAR})
	file(MD5 ${MODELS_DIR}/${CLASS_RANKER_MODEL_FILE_NAME_TAR} CHECKSUM)
	if (${CHECKSUM} STREQUAL ${CLASS_RANKER_MODEL_MD5})
		set(CLASS_RANKER_MODEL_FOUND 1)
		message(STATUS "Found local copy of class ranker model")
	else()
		message(STATUS "Checksum of local class ranker model did not match...")
	endif()
endif()

if(NOT CLASS_RANKER_MODEL_FOUND)
	externalproject_add(class_ranker_model_file
			URL ${CLASS_RANKER_MODEL_URL}
			URL_MD5 ${CLASS_RANKER_MODEL_MD5}
			DOWNLOAD_DIR ${MODELS_DIR}
			SOURCE_DIR "${MODELS_DIR}/class_ranker"
			BUILD_COMMAND ""
			INSTALL_COMMAND ""
			CONFIGURE_COMMAND ""
			#	LOG_CONFIGURE
			#	LOG_BUILD
			LOG_INSTALL)

	message(STATUS "--------------------------------------------------------")
	message(STATUS "--------- FOUND NO LOCAL COPY OF TORCH MODELS. ---------")
	message(STATUS "------- WILL BE DOWNLOADED DURING COMPILE-TIME. --------")
	message(STATUS "--------------------------------------------------------")
	message(STATUS "---- A WORKING INTERNET CONNECTION WILL BE REQUIRED. ---")
	message(STATUS "-- TO SKIP, RECONFIGURE WITH -DFETCH_TORCH_MODELS=OFF --")
	message(STATUS "--------------------------------------------------------")
else()
	add_custom_target(class_ranker_model_file ALL)
endif()

add_custom_command(TARGET class_ranker_model_file POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy
		${MODELS_DIR}/class_ranker/${CLASS_RANKER_MODEL_FILE_NAME}
		${CMAKE_BINARY_DIR}/bin/${CLASS_RANKER_DEFAULT_MODEL_FILE_NAME})

#add_definitions(-DCLASS_RANKER_DEFAULT_MODEL_FILE_NAME="${CLASS_RANKER_DEFAULT_MODEL_FILE_NAME}")
