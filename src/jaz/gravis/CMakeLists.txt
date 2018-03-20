add_subdirectory(colour)
add_subdirectory(private)
add_subdirectory(tImage)
add_subdirectory(io)

set( install_files 
    Exception.h
    lapack.h
    matrix_blas_reference.h
    Mesh.h
    NMesh.h
    NTuple.h
    OBJReader.h
    OBJWriter.h
    program_options.h
    StringFormat.h
    t2Matrix.h
    t2Vector.h
    t3Matrix.h
    t3Vector.h
    t4Matrix.h
    t4Vector.h
    tArray.h
    tBGR.h
    tDefaultVector.h
    tGray_A.h
    tImageAlgorithm.h
    tImage.h
    Timer.hpp
    tLab.h
    tMatrix.h
    tMesh.h
    tMM.h
    tQuaternion.h
    tRGB_A.h
    tRGBA.h
    tRGB.h
    Tuple.h
    tVarMatrix_blas.h
    tVarMatrix_blas.hxx
    tVarMatrix.h
    tYCbCr.h
    PushPull.h
)
INSTALL(FILES ${install_files} 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${LIBTITLE}-${LIBVERSION}/${LIBTITLE})
