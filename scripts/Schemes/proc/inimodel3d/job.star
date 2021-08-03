
# version 30001

data_job

_rlnJobTypeLabel             relion.initialmodel
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
ctf_intact_first_peak         No 
do_combine_thru_disc         No 
do_ctf_correction        Yes 
do_parallel_discio        Yes 
do_preread_images         No 
  do_queue         No 
do_solvent        Yes 
   fn_cont         "" 
    fn_img Schemes/proc/select_parts/particles.star 
   gpu_ids        0,1 
min_dedicated         24 
nr_classes          3 
   nr_iter        200 
    nr_mpi          1 
   nr_pool         30 
nr_threads         12 
offset_range          6 
offset_step          2 
other_args         "" 
particle_diameter        200 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
  sampling "15 degrees" 
scratch_dir       "" 
skip_gridding        Yes 
  sym_name         C1
 tau_fudge          4 
   use_gpu        Yes 
 
