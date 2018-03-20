set( install_files
    absolute_paths.h
    msh.h
    obj.h
    OBJReader.h
    OBJWriter.h
    )
INSTALL(FILES ${install_files} 
        DESTINATION
        ${CMAKE_INSTALL_PREFIX}/include/${LIBTITLE}-${LIBVERSION}/${LIBTITLE}/io/mesh)
