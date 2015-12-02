
if(FORCE_OWN_FLTK AND ALLOW_OWN_FLTK)       
    set(FLTK_EXTERNAL_PATH "${CMAKE_SOURCE_DIR}/external/fltk")
    set(ext_conf_flags_fltk --enable-shared --prefix=${FLTK_EXTERNAL_PATH})
    
    ## ------------------------------------------------------------- PREVIOUS EXT LIBS? --
    
    find_path(FLTK_PATH         NAMES fl.h         PATHS ${FLTK_EXTERNAL_PATH}/include NO_DEFAULT_PATH)
    find_path(FLTK_INCLUDES     NAMES fl.h         PATHS ${FLTK_EXTERNAL_PATH}/include NO_DEFAULT_PATH) 
    find_library(FLTK_LIBRARIES NAMES libfltk.so   PATHS ${FLTK_EXTERNAL_PATH}/lib NO_DEFAULT_PATH)   
    
    if(FLTK_PATH AND FLTK_INCLUDES AND FLTK_LIBRARIES)
        set(FLTK_FOUND TRUE)
        message( STATUS "found previously built external (non-system) fltk lib")
    endif()
    ## ----------------------------------------------------------------- NEW EXT LIBS? --  
     
    if(NOT FLTK_FOUND)
    
        include(ExternalProject)
        set(FLTK_EXTERNAL_LIBS_TAR_DIRECTORY  ${FLTK_EXTERNAL_PATH})
        set(FLTK_EXTERNAL_LIBS_EXTRACT_TARGET ${FLTK_EXTERNAL_LIBS_TAR_DIRECTORY})
        message( STATUS "no fltk found, the following paths are set for libs/headers TO BE built")
        
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
#       TIMEOUT 15
        URL_MD5 ${FLTK_HASH}
        DOWNLOAD_DIR ${FLTK_EXTERNAL_LIBS_TAR_DIRECTORY}
        DOWNLOAD_NAME ${FLTK_TAR_NAME}
        SOURCE_DIR ${FLTK_LIB_DIR}
        CONFIGURE_COMMAND <SOURCE_DIR>/configure ${ext_conf_flags_fltk}
        INSTALL_DIR ${FLTK_EXTERNAL_PATH}/fltk
        BINARY_DIR ${FLTK_EXTERNAL_PATH}/fltk
        BUILD_COMMAND ${MAKE}
#	    LOG_CONFIGURE
#	    LOG_BUILD
        LOG_INSTALL)
        
        set(NEW_OWN_FLTK TRUE)
        set(FLTK_FOUND TRUE)
        
        set(FLTK_LIBRARIES     "${FLTK_EXTERNAL_PATH}/lib/libfltk.so" )  
#       set(FLTK_PATH          "${FLTK_EXTERNAL_PATH}/includes/FL/Fl.H" )  
#       set(FLTK_INCLUDES      "${FLTK_EXTERNAL_PATH}/includes/FL/Fl.H" ) 
        set(FLTK_INCLUDE_DIR   "${FLTK_EXTERNAL_PATH}/includes" )
    
    endif()
else()
    message( STATUS "\n-- ------------------ YOU HAVE NO FLTK-LIBS ------------------")
    message( STATUS "CCmake found no fltk-libs on your system, which are required for the GUI. Either ")
    message( STATUS "     a) Make sure cmake can find them ")
    message( STATUS "     b) Add the flag -DOWN_FLTK=ON to your cmake command to build a local fltk-lib")
    message( STATUS "     c) Add the flag -DGUI=OFF to avoid using fltk" )
    message( STATUS "--------------------------------------------------------")
	message( FATAL_ERROR "FLTK is required." )
endif()
