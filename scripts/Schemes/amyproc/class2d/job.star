
# version 50001

data_job

_rlnJobTypeLabel             relion.class2d.helical
_rlnJobIsContinue                       0
_rlnJobIsTomo                           0
 

# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
allow_coarser         No 
ctf_intact_first_peak        Yes 
do_bimodal_psi        Yes 
 do_center         No 
do_combine_thru_disc         No 
do_ctf_correction        Yes 
     do_em        Yes 
   do_grad         No 
  do_helix        Yes 
do_parallel_discio        Yes 
do_preread_images         No 
  do_queue         No 
do_restrict_xoff        Yes 
do_zero_mask        Yes 
dont_skip_align        Yes 
   fn_cont         "" 
    fn_img $$mybatch
   gpu_ids    0:1:2:3 
helical_rise       4.75 
helical_tube_outer_diameter        200 
highres_limit         -1 
min_dedicated         24 
nr_classes        100 
nr_iter_em         25 
nr_iter_grad        200 
    nr_mpi          5 
   nr_pool         30 
nr_threads          6 
offset_range          6 
offset_step          1 
other_args         "" 
particle_diameter        600 
psi_sampling          2 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-gpu-4.0.csh 
 queuename    openmpi 
 range_psi          6 
scratch_dir       /ssd 
 tau_fudge          2 
   use_gpu        Yes 
 
