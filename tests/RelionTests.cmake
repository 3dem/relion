file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/test_output)

#--------------------   CPU-computation test  -----------------------
#--------------------------------------------------------------------
add_test(NAME deleting_old_files
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
        COMMAND rm gpu_out_10k_diff2s.txt cpu_out_10k_diff2s.txt )
        
add_test(NAME CPU-3Dc_produce_data
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 40
        --ctf 
        --ctf_corrected_ref
        --iter 2 
        --tau2_fudge 2 
        --K 1 
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
        --onthefly_shifts
        --scale)     
SET_TESTS_PROPERTIES(CPU-3Dc_produce_data PROPERTIES DEPENDS deleting_old_files)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        cpu_out_shifted_image.mrc)
        SET_TESTS_PROPERTIES(CPU-3Dc_shifted_image PROPERTIES DEPENDS CPU-3Dc_produce_data_and_diff2)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_10kdiffs
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_diff2s.txt 
        cpu_out_10k_diff2s.txt)
SET_TESTS_PROPERTIES(CPU-3Dc_10kdiffs PROPERTIES DEPENDS CPU-3Dc_produce_data)
#--------------------------------------------------------------------  
#add_test(NAME CPU-3Dc_ref
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#         ref_3D_fref.mrc 
#        cpu_out_fref.mrc)
#        SET_TESTS_PROPERTIES(CPU-3Dc_ref PROPERTIES DEPENDS CPU-3Dc_produce_data_and_diff2)           
#--------------------------------------------------------------------
#add_test(NAME CPU-3Dc_ctf
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#        ref_3D_ctf.mrc 
#       cpu_out_ctf.mrc)
#       SET_TESTS_PROPERTIES(CPU-3Dc_ctf PROPERTIES DEPENDS CPU-3Dc_produce_data_and_diff2)
#--------------------------------------------------------------------
#add_test(NAME CPU-3Dc_refctf
#       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#       COMMAND ${CMAKE_COMMAND} -E compare_files 
#        ref_3D_frefctf.mrc 
#       cpu_out_frefctf.mrc)
#       SET_TESTS_PROPERTIES(CPU-3Dc_refctf PROPERTIES DEPENDS CPU-3Dc_produce_data_and_diff2)
#--------------------------------------------------------------------
#add_test(NAME CPU-3Dc_ref-vs-refctf
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#        cpu_out_ref.mrc 
#        cpu_out_frefctf.mrc)
#SET_TESTS_PROPERTIES(CPU-3Dc_ref-vs-refctf PROPERTIES WILL_FAIL TRUE DEPENDS CPU-3Dc_produce_data_and_diff2)


#--------------------   GPU-computation test  -----------------------
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_K1-produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 40 
        --ctf
        --ctf_corrected_ref
        --iter 2 
        --tau2_fudge 2 
        --K 1 
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
        --onthefly_shifts
        --gpu
        --scale)
#-------------------------------------------------------------------- 
add_test(NAME GPU-3Dc_K1-shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        gpu_out_shifted_image.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_K1-shifted_image PROPERTIES DEPENDS GPU-3Dc_K1-produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_K1-10kdiffs
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_diff2s.txt 
        gpu_out_10k_diff2s.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_K1-10kdiffs PROPERTIES DEPENDS GPU-3Dc_K1-produce)        
#--------------------------------------------------------------------


#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_K4-produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 40 
        --ctf
        --ctf_corrected_ref
        --iter 2 
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
        --onthefly_shifts
        --gpu)
SET_TESTS_PROPERTIES(GPU-3Dc_K4-produce PROPERTIES DEPENDS GPU-3Dc_K1-10kdiffs) 
#-------------------------------------------------------------------- 
add_test(NAME GPU-3Dc_K4-shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        gpu_out_shifted_image.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_K4-shifted_image PROPERTIES DEPENDS GPU-3Dc_K4-produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_K4-10kdiffs
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_diff2s.txt 
        gpu_out_10k_diff2s.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_K4-10kdiffs PROPERTIES DEPENDS GPU-3Dc_K4-produce)           
#-------------------------------------------------------------------- 

add_test(NAME GPU-3Dc_noCTF-produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 40 
        --ctf_corrected_ref
        --iter 2 
        --tau2_fudge 2 
        --K 1 
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
        --onthefly_shifts
        --gpu)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF-produce PROPERTIES DEPENDS GPU-3Dc_K4-10kdiffs) 
#-------------------------------------------------------------------- 
add_test(NAME GPU-3Dc_noCTF-shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        gpu_out_shifted_image.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF-shifted_image PROPERTIES DEPENDS GPU-3Dc_noCTF-produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF-10kdiffs
WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_diff2s_noCTF.txt 
        gpu_out_10k_diff2s.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF-10kdiffs PROPERTIES DEPENDS GPU-3Dc_noCTF-produce)        
#--------------------------------------------------------------------
           
