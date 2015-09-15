file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/test_output)

add_test(NAME deleting_old_files
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests
        COMMAND rm gpu_out_10k_weights.txt 
				   cpu_out_10k_weigths.txt
				   gpu_out_10k_weights_afterstore.txt 
				   cpu_out_10k_weigths_afterstore.txt )
				   
set(RELION_INPUT_DIR_RYR "/nethome/projects/RELION/refData/inputData/ryr")
			   
#--------------------------------------------------------------------  
add_test(NAME CPU-3Dc_tut
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu3Dc_tut
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
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_tut
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_tut
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
add_test(NAME GPU-3Dc_noCTF_tut
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_NoCTF_tut
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
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_10P-ryr_K1
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_ryr
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_10_local.star" 
        --j 1 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 8 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 1 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 3 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --memory_per_thread 2 
        --random_seed 1 
        --onthefly_shifts 
        --gpu) 
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_2kP-ryr_K3
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_ryr
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_2k_local.star" 
        --j 1 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 8 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 3 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 3 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --memory_per_thread 2 
        --random_seed 1 
        --onthefly_shifts 
        --gpu) 
#--------------------------------------------------------------------
add_test(NAME GPU-2Dc_100P-empiar_K10
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu2Dc_s_empiar_100 
        --i /nethome/projects/RELION/refData/inputData/empiar/particles_rh15_100_L.star 
        --particle_diameter 200 
        --angpix 1.77 
        --iter 3 
        --tau2_fudge 2 
        --K 10 
        --ctf 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --psi_step 10 
        --offset_range 5 
        --offset_step 2 
        --norm 
        --scale 
        --j 1 
        --onthefly_shifts 
        --memory_per_thread 8 
        --random_seed 1993 
        --gpu)
