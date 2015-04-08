file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/test_output)

#--------------------   CPU-computation test  -----------------------
#--------------------------------------------------------------------
add_test(NAME deleting_old_files
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
        COMMAND rm gpu_out_10k_weights.txt 
				   cpu_out_10k_weigths.txt
				   gpu_out_10k_weights_afterstore.txt 
				   cpu_out_10k_weigths_afterstore.txt )
        
add_test(NAME CPU-3Dc_produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 50
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
SET_TESTS_PROPERTIES(CPU-3Dc_produce PROPERTIES DEPENDS deleting_old_files)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        cpu_out_shifted_image.mrc)
        SET_TESTS_PROPERTIES(CPU-3Dc_shifted_image PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_ctf
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_ctf.mrc 
       cpu_out_ctf.mrc)
       SET_TESTS_PROPERTIES(CPU-3Dc_ctf PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------


#---------WILL IMPLEMENT THESE SOON---------------------
#add_test(NAME CPU-3Dc_ref
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#         ref_3D_fref.mrc 
#        cpu_out_fref.mrc)
#        SET_TESTS_PROPERTIES(CPU-3Dc_ref PROPERTIES DEPENDS CPU-3Dc_produce)           
#--------------------------------------------------------------------
#add_test(NAME CPU-3Dc_refctf
#       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#       COMMAND ${CMAKE_COMMAND} -E compare_files 
#        ref_3D_frefctf.mrc 
#       cpu_out_frefctf.mrc)
#       SET_TESTS_PROPERTIES(CPU-3Dc_refctf PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------


add_test(NAME CPU-3Dc_10kweights
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_weights.txt 
        cpu_out_10k_weights.txt)
SET_TESTS_PROPERTIES(CPU-3Dc_10kweights PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_10kweights_afterstore
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_weights_afterstore.txt 
        cpu_out_10k_weights_afterstore.txt)
SET_TESTS_PROPERTIES(CPU-3Dc_10kweights_afterstore PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_mresol_coarse
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_coarse.mrc
        cpu_out_mresol_coarse.mrc)
SET_TESTS_PROPERTIES(CPU-3Dc_mresol_coarse PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_mresol_fine
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_fine.mrc
        cpu_out_mresol_fine.mrc)
SET_TESTS_PROPERTIES(CPU-3Dc_mresol_fine PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_norm_correction
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_norm_correction.mrc
        cpu_out_norm_correction.mrc)
SET_TESTS_PROPERTIES(CPU-3Dc_norm_correction PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_sigma2_noise
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_sigma2_noise.mrc
        cpu_out_sigma2_noise.mrc)
SET_TESTS_PROPERTIES(CPU-3Dc_sigma2_noise PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------


  
#--------------------------------------------------------------------
#--------------------   GPU-computation test  -----------------------
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc-produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 50 
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
add_test(NAME GPU-3Dc-shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_shifted_image.mrc 
        gpu_out_shifted_image.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc-shifted_image PROPERTIES DEPENDS GPU-3Dc-produce)      
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_ctf
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_3D_ctf.mrc 
       gpu_out_ctf.mrc)
       SET_TESTS_PROPERTIES(GPU-3Dc_ctf PROPERTIES DEPENDS GPU-3Dc-produce)
#--------------------------------------------------------------------

#---------WILL IMPLEMENT THESE SOON---------------------
#add_test(NAME GPU-3Dc_ref
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#         ref_3D_fref.mrc 
#        gpu_out_fref.mrc)
#        SET_TESTS_PROPERTIES(GPU-3Dc_ref PROPERTIES DEPENDS GPU-3Dc-produce)           

#--------------------------------------------------------------------
#add_test(NAME GPU-3Dc_refctf
#       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#       COMMAND ${CMAKE_COMMAND} -E compare_files 
#        ref_3D_frefctf.mrc 
#       gpu_out_frefctf.mrc)
#       SET_TESTS_PROPERTIES(GPU-3Dc_refctf PROPERTIES DEPENDS GPU-3Dc-produce)


add_test(NAME GPU-3Dc_10kweights
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_weights.txt 
        gpu_out_10k_weights.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_10kweights PROPERTIES DEPENDS GPU-3Dc-produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_10kweights_afterstore
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_10k_weights_afterstore.txt 
        gpu_out_10k_weights_afterstore.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_10kweights_afterstore PROPERTIES DEPENDS GPU-3Dc-produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_mresol_coarse
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_coarse.mrc
        gpu_out_mresol_coarse.mrc)
SET_TESTS_PROPERTIES(CPU-3Dc_mresol_coarse PROPERTIES DEPENDS CPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_mresol_fine
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_fine.mrc
        gpu_out_mresol_fine.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_mresol_fine PROPERTIES DEPENDS GPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_norm_correction
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_norm_correction.mrc
        cpu_out_norm_correction.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_norm_correction PROPERTIES DEPENDS GPU-3Dc_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_sigma2_noise
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_sigma2_noise.mrc
        cpu_out_sigma2_noise.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_sigma2_noise PROPERTIES DEPENDS GPU-3Dc_produce)
#-------------------------------------------------------------------- 

#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF-produce
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/out
        --i ap_s_c2d_red.star 
        --particle_diameter 200 
        --angpix 1.0 
        --ref ref_model.mrc 
        --ini_high 50 
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
        SET_TESTS_PROPERTIES(GPU-3Dc_noCTF-produce PROPERTIES DEPENDS GPU-3Dc_mresol_fine)  
#-------------------------------------------------------------------- 
add_test(NAME GPU-3Dc_noCTF-shifted_image
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_noCTF_shifted_image.mrc 
          gpu_out_shifted_image.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF-shifted_image PROPERTIES DEPENDS GPU-3Dc_noCTF-produce)      
#--------------------------------------------------------------------#--------------------------------------------------------------------

#---------WILL IMPLEMENT THESE SOON---------------------
#add_test(NAME GPU-3Dc_noCTF_ref
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
#        COMMAND ${CMAKE_COMMAND} -E compare_files 
#        ref_noCTF_3D_fref.mrc 
#             gpu_out_fref.mrc)
#        SET_TESTS_PROPERTIES(GPU-3Dc_noCTF_ref PROPERTIES DEPENDS GPU-3Dc_noCTF-produce)  
        
                 
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF_10kweights
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
      ref_noCTF_10k_weights.txt 
        gpu_out_10k_weights.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF_10kweights PROPERTIES DEPENDS GPU-3Dc_noCTF_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF_10kweights_afterstore
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
        ref_noCTF_10k_weights_afterstore.txt 
          gpu_out_10k_weights_afterstore.txt)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF_10kweights_afterstore PROPERTIES DEPENDS GPU-3Dc_noCTF_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF_mresol_coarse
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_coarse.mrc
        gpu_out_mresol_coarse.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF_mresol_coarse PROPERTIES DEPENDS GPU-3Dc_noCTF_produce)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_noCTF_mresol_fine
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests
        COMMAND ${CMAKE_COMMAND} -E compare_files 
            ref_mresol_fine.mrc
        gpu_out_mresol_fine.mrc)
SET_TESTS_PROPERTIES(GPU-3Dc_noCTF_mresol_fine PROPERTIES DEPENDS GPU-3Dc_noCTF_produce)
#--------------------------------------------------------------------
