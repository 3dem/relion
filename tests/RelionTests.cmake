file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/test_output)
#--------------------------------------------------------------------
add_test(NAME 3Dclassification_produce_data_and_diff2 
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1 
        --ref ref_model.mrc 
        --ini_high 50 
        --ctf 
        --ctf_corrected_ref
        --iter 1 
        --tau2_fudge 2 
        --K 4 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 2 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --j 1 
        --memory_per_thread 4 
        --random_seed 1993 
        --onthefly_shifts)      
#--------------------------------------------------------------------
add_test(NAME 3Dc_ref
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_fref.mrc 
           out_fref.mrc)           
#--------------------------------------------------------------------
add_test(NAME 3Dc_ctf
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_ctf.mrc 
           out_ctf.mrc)
#--------------------------------------------------------------------
add_test(NAME 3Dc_refctf
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_frefctf.mrc 
           out_frefctf.mrc)
#--------------------------------------------------------------------
add_test(NAME 3Dc_ref-vs-refctf
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
           out_ref.mrc 
           out_frefctf.mrc)
SET_TESTS_PROPERTIES(3Dc_ref-vs-refctf PROPERTIES WILL_FAIL TRUE)
#--------------------------------------------------------------------
add_test(NAME 3Dc_otfshift
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_otfshift.mrc 
           out_otfshift.mrc)
           
