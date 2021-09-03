
# version 30001

data_job

_rlnJobTypeLabel             relion.refine3d
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
auto_faster        No
auto_local_sampling "1.8 degrees" 
ctf_intact_first_peak         No 
do_apply_helical_symmetry        Yes 
do_combine_thru_disc         No 
do_ctf_correction        Yes 
  do_helix         No 
do_local_search_helical_symmetry         No 
   do_pad1         No 
do_parallel_discio        Yes 
do_preread_images         No 
  do_queue         No 
do_solvent_fsc         No 
do_zero_mask        Yes 
   fn_cont         "" 
    fn_img Schemes/proc/select_parts/particles.star
   fn_mask         "" 
    fn_ref $$myref
   gpu_ids        0:1 
helical_nr_asu          1 
helical_range_distance         -1 
helical_rise_inistep          0 
helical_rise_initial          0 
helical_rise_max          0 
helical_rise_min          0 
helical_tube_inner_diameter         -1 
helical_tube_outer_diameter         -1 
helical_twist_inistep          0 
helical_twist_initial          0 
helical_twist_max          0 
helical_twist_min          0 
helical_z_percentage         30 
  ini_high         40 
keep_tilt_prior_fixed        Yes 
min_dedicated         24 
    nr_mpi          3 
   nr_pool         30 
nr_threads         8 
offset_range          5 
offset_step          1 
other_args         "" 
particle_diameter        200 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
 range_psi         10 
 range_rot         -1 
range_tilt         15 
ref_correct_greyscale        Yes 
 relax_sym         "" 
  sampling "7.5 degrees" 
scratch_dir       "" 
skip_gridding        Yes 
  sym_name         C1
   use_gpu        Yes 
 
