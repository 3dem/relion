set( install_files
    rgb_hsv.h
    )
INSTALL(FILES ${install_files} 
        DESTINATION
        ${CMAKE_INSTALL_PREFIX}/include/${LIBTITLE}-${LIBVERSION}/${LIBTITLE}/colour)
