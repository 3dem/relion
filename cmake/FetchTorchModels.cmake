#################################################################################
# Description:  Cmake helper to download and configure LibTorch models.
# Author:       Dari Kimanius (github.com/dkimanius)
# Date created: 2021-01-14
#################################################################################

include(ExternalProject)
set(CLASS_RANKER_DEFAULT_MODEL_FILE_NAME "relion_class_ranker_default_model.pt")

set(MODELS_DIR  ${CMAKE_SOURCE_DIR}/external/torch_models)
if(NOT EXISTS ${MODELS_DIR})
	file(MAKE_DIRECTORY ${MODELS_DIR})
endif()

set(CLASS_RANKER_MODEL_FILE_NAME class_ranker_0.1_batch32_rate5e-05_drop0.3_epoch50_model.pt)
set(CLASS_RANKER_MODEL_GZIP_FILE_NAME ${CLASS_RANKER_MODEL_FILE_NAME}.gz)
set(CLASS_RANKER_MODEL_URL "ftp://ftp.mrc-lmb.cam.ac.uk/pub/dari/${CLASS_RANKER_MODEL_GZIP_FILE_NAME}")
set(CLASS_RANKER_MODEL_MD5 f699280eb98bcbcb693aea9a266c1a23)

message(STATUS "CLASS_RANKER_MODEL_GZIP_FILE_NAME ${CLASS_RANKER_MODEL_GZIP_FILE_NAME}")

externalproject_add(class_ranker_model_file
		URL ${CLASS_RANKER_MODEL_URL}
		URL_MD5 ${CLASS_RANKER_MODEL_MD5}
		DOWNLOAD_DIR ${MODELS_DIR}
		BUILD_COMMAND ""
		INSTALL_COMMAND ""
		DOWNLOAD_NO_EXTRACT 1 # Pytorch PT files are tar archives
		DOWNLOAD_NO_PROGRESS 1
		CONFIGURE_COMMAND gzip -fd <DOWNLOADED_FILE>
		#	LOG_CONFIGURE
		#	LOG_BUILD
		LOG_INSTALL)

add_custom_command(TARGET class_ranker_model_file POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy
		${MODELS_DIR}/${CLASS_RANKER_MODEL_FILE_NAME}
		${CMAKE_BINARY_DIR}/bin/${CLASS_RANKER_DEFAULT_MODEL_FILE_NAME})

add_definitions(-DCLASS_RANKER_DEFAULT_MODEL_FILE_NAME="${CLASS_RANKER_DEFAULT_MODEL_FILE_NAME}")