
# version 50001

data_job

_rlnJobTypeLabel             relion.motioncorr
_rlnJobIsContinue                       0
_rlnJobIsTomo                           0
 

# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
   bfactor        150 
bin_factor          1 
do_dose_weighting        Yes 
do_float16        Yes 
do_own_motioncor        Yes 
  do_queue         No 
do_save_noDW         No 
do_save_ps        Yes 
dose_per_frame        1.0 
eer_grouping         -1 
first_frame_sum          1 
 fn_defect         "" 
fn_gain_ref   gain.mrc 
fn_motioncor2_exe /public/EM/MOTIONCOR2/MotionCor2 
 gain_flip "No flipping (0)" 
  gain_rot "No rotation (0)" 
   gpu_ids          0 
group_for_ps          4 
group_frames          1 
input_star_mics Schemes/amyprep/import/movies.star 
last_frame_sum         -1 
min_dedicated         24 
    nr_mpi         9 
nr_threads         12 
other_args "--do_at_most $$do_at_most" 
other_motioncor2_args         "" 
   patch_x          4 
   patch_y          4 
pre_exposure          0 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-cpu-devel.sh 
 queuename    openmpi 
 
