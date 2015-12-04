message(STATUS "-------------------------------------------------")   
message(STATUS "------------- USING LOCAL FLTK LIBS -------------")  
message(STATUS "-------------------------------------------------") 

set(FLTK_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fltk")
set(ext_conf_flags_fltk --enable-shared --prefix=${FLTK_EXTERNAL_PATH})

## ------------------------------------------------------------- PREVIOUS EXT LIBS? --

find_library(FLTK_LIBRARIES NAMES libfltk.so  PATHS  "${FLTK_EXTERNAL_PATH}/lib" NO_DEFAULT_PATH) 
find_path(FLTK_INCLUDE_DIR  NAMES Fl.H   PATHS  "${FLTK_EXTERNAL_PATH}/include/FL" NO_DEFAULT_PATH) 
find_path(FLTK_PATH         NAMES Fl.H   PATHS  "${FLTK_EXTERNAL_PATH}/include/FL" NO_DEFAULT_PATH) 
find_path(FLTK_INCLUDES     NAMES Fl.H   PATHS  "${FLTK_EXTERNAL_PATH}/include/FL" NO_DEFAULT_PATH)  

if(FLTK_INCLUDE_DIR AND FLTK_LIBRARIES)
    set(FLTK_FOUND TRUE)
    message( STATUS "found previously built external (non-system) fltk lib")
endif()
## ----------------------------------------------------------------- NEW EXT LIBS? --  
 
if(NOT FLTK_FOUND)

    set(FLTK_LIBRARIES     "${FLTK_EXTERNAL_PATH}/lib/libfltk.so" ) 
    set(FLTK_INCLUDE_DIR   "${FLTK_EXTERNAL_PATH}/includes" )
    #set(FLTK_PATH          "${FLTK_EXTERNAL_PATH}/includes/FL/Fl.H" )  
    #set(FLTK_INCLUDES      "${FLTK_EXTERNAL_PATH}/includes/FL/Fl.H" ) 

    include(ExternalProject)
    set(FLTK_EXTERNAL_LIBS_TAR_DIRECTORY  ${FLTK_EXTERNAL_PATH})
    set(FLTK_EXTERNAL_LIBS_EXTRACT_TARGET ${FLTK_EXTERNAL_LIBS_TAR_DIRECTORY})
       
    message(STATUS "no previous fltk found, the following paths are set for libs/headers TO BE built")
    
    #set(FLTK_TAR_FILE https://drive.google.com/uc?export=download&id=0B942d76zVnSeUWgyaklWOFZlN2s)  # FLTK 1.3.0
    #set(FLTK_HASH cf3687ed404bd72347466663b90a8b09)
    #set(FLTK_TAR_NAME fltk-1.3.0.tar.gz)

    set(FLTK_TAR_FILE https://drive.google.com/uc?export=download&id=0B942d76zVnSeazZWcExRaXIyVDg)   # FLTK 1.3.3
    set(FLTK_HASH 9ccdb0d19dc104b87179bd9fd10822e3)
    set(FLTK_TAR_NAME fltk-1.3.3-source.tar.gz)
    
    set(FLTK_LIB_DIR ${FLTK_EXTERNAL_LIBS_EXTRACT_TARGET}/fltk)
    set(FLTK_BUILD_DIR ${FLTK_EXTERNAL_LIBS_EXTRACT_TARGET}/fltk-build)
   
    set(CMAKE_INSTALL_PREFIX  ${FLTK_EXTERNAL_PATH})
    externalproject_add(FLTK
    URL ${FLTK_TAR_FILE}
#   TIMEOUT 15
    URL_MD5 ${FLTK_HASH}
    DOWNLOAD_DIR ${FLTK_EXTERNAL_LIBS_TAR_DIRECTORY}
    DOWNLOAD_NAME ${FLTK_TAR_NAME}
    SOURCE_DIR ${FLTK_LIB_DIR}
    CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fltk}
    INSTALL_DIR ${FLTK_EXTERNAL_PATH}/fltk
    BINARY_DIR ${FLTK_EXTERNAL_PATH}/fltk
    BUILD_COMMAND ${MAKE}
#	LOG_CONFIGURE
#	LOG_BUILD
    LOG_INSTALL)
    
    set(NEW_OWN_FLTK TRUE)
    set(FLTK_FOUND TRUE)

endif()

message(STATUS "FLTK_INCLUDE_DIR: ${FLTK_INCLUDE_DIR}")
message(STATUS "FLTK_LIBRARIES:   ${FLTK_LIBRARIES}")

