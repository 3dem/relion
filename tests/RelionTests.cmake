file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/test_output)

set(RELION_INPUT_DIR_RYR "/nethome/projects/RELION/refData/inputData/ryr")
set(RELION_INPUT_DIR_EMP "/nethome/projects/RELION/refData/inputData/empiar")
set(RELION_INPUT_DIR_TAG "/nethome/projects/RELION/refData/inputData/tagcase")

#--------------------------------------------------------------------			   
add_test(NAME GPU-2Dc_small
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu2Dc_small
        --i "${RELION_INPUT_DIR_EMP}/particles_rh15_100_L.star" 
        --particle_diameter 200 
        --angpix 1.77 
        --iter 3 
        --tau2_fudge 2 
        --K 10 
        --ctf 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --psi_step 11.25
        --offset_range 5 
        --offset_step 2 
        --norm 
        --scale 
        --j 1 
        --onthefly_shifts 
        --random_seed 1993 
	--perturb 0
        --gpu)
#--------------------------------------------------------------------
add_test(NAME CPU-2Dc_small
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu2Dc_small
        --i "${RELION_INPUT_DIR_EMP}/particles_rh15_100_L.star" 
        --particle_diameter 200 
        --angpix 1.77 
        --iter 3 
        --tau2_fudge 2 
        --K 10 
        --ctf 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --psi_step 11.25
        --offset_range 5 
        --offset_step 2 
        --norm 
        --scale 
        --j 1 
        --onthefly_shifts 
        --random_seed 1993 
	--perturb 0)
#--------------------------------------------------------------------
add_test(NAME GPU-2Dc_medium
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu2Dc_medium 
        --i "${RELION_INPUT_DIR_EMP}/particles_rh15_600_L.star" 
        --particle_diameter 200 
        --angpix 1.77 
        --iter 5 
        --tau2_fudge 2 
        --K 50 
        --ctf 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --psi_step 5.625
        --offset_range 5 
        --offset_step 1 
        --norm 
        --scale 
        --j 1 
        --random_seed 1993 
	--perturb 0
        --gpu)
#--------------------------------------------------------------------
add_test(NAME CPU-2Dc_medium
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu2Dc_medium 
        --i "${RELION_INPUT_DIR_EMP}/particles_rh15_600_L.star" 
        --particle_diameter 200 
        --angpix 1.77 
        --iter 5 
        --tau2_fudge 2 
        --K 50 
        --ctf 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --psi_step 5.625
        --offset_range 5 
        --offset_step 1 
        --norm 
        --scale 
        --j 5 
        --random_seed 1993 
	--perturb 0)
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_small
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_small
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_10_local.star" 
        --j 1 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 1 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 2 
        --offset_range 3 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1 
        --onthefly_shifts 
	--perturb 0
        --gpu) 
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_small
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu3Dc_small
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_10_local.star" 
        --j 1 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 1 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 2 
        --offset_range 3 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1 
        --onthefly_shifts 
	--perturb 0) 
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_medium
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_medium
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_500_local.star" 
        --j 5 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 3 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 2 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1
        --onthefly_shifts 
	--perturb 0
        --gpu) 
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_medium
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu3Dc_medium
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_500_local.star" 
        --j 5 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 3 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 2 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1
        --onthefly_shifts 
	--perturb 0) 
#--------------------------------------------------------------------
add_test(NAME GPU-3Dc_large
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu3Dc_large
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_2k_local.star" 
        --j 5 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 5 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 3 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1 
        --onthefly_shifts 
	--perturb 0
        --gpu) 
#--------------------------------------------------------------------
add_test(NAME CPU-3Dc_large
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu3Dc_large
        --i "${RELION_INPUT_DIR_RYR}/benchmark_ryr_2k_local.star" 
        --j 5 
        --particle_diameter 400 
        --angpix 1.34 
        --ref "${RELION_INPUT_DIR_RYR}/ref_8A.mrc" 
        --ini_high 30 
        --ctf 
        --ctf_corrected_ref 
        --iter 2 
        --tau2_fudge 2 
        --K 5 
        --flatten_solvent 
        --zero_mask 
        --oversampling 1 
        --healpix_order 3 
        --offset_range 5 
        --offset_step 2 
        --sym C1 
        --norm 
        --scale 
        --random_seed 1 
        --onthefly_shifts 
	--perturb 0) 
#--------------------------------------------------------------------
add_test(NAME CPU-2Dc_tag_helix
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/cpu_2Dc_tag_helix
        --i "${RELION_INPUT_DIR_TAG}/helixtest.star" 
        --j 5 
        --particle_diameter 160 
        --angpix 2.8 
        --iter 10 
        --tau2_fudge 4 
        --K 2 
        --flatten_solvent 
        --psi_step 5.625 
        --zero_mask 
        --oversampling 1 
        --offset_range 5 
        --offset_step 2 
        --scale 
        --random_seed 1 
        --helical_outer_diameter 100
        --bimodal_psi
        --sigma_psi 3.3333) 
        #--------------------------------------------------------------------
add_test(NAME GPU-2Dc_tag_helix
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/tests 
        COMMAND refine --o test_output/gpu_2Dc_tag_helix
        --i "${RELION_INPUT_DIR_TAG}/helixtest.star" 
        --j 5 
        --particle_diameter 160 
        --angpix 2.8 
        --iter 10 
        --tau2_fudge 4 
        --K 2 
        --flatten_solvent 
        --psi_step 5.625 
        --zero_mask 
        --oversampling 1 
        --offset_range 5 
        --offset_step 2 
        --scale 
        --random_seed 1 
        --helical_outer_diameter 100
        --bimodal_psi
        --sigma_psi 3.3333
        --gpu) 
	
