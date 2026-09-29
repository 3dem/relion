
# version 50001

data_job

_rlnJobTypeLabel             relion.autopick.amypick
_rlnJobIsContinue                       1
_rlnJobIsTomo                           0
 

# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
amyloid_threshold        .25 
    angpix         -1 
angpix_ref         -1 
continue_manual         No 
do_amyloid        Yes 
do_amyloid_fom        Yes 
do_amyloid_plot         No 
do_amyloid_tracing         No 
do_ctf_autopick        Yes 
do_ignore_first_ctfpeak_autopick         No 
do_invert_refs        Yes 
    do_log         No 
do_pick_helical_segments        Yes 
  do_queue         No 
do_read_fom_maps         No 
  do_ref3d         No 
   do_refs         No 
  do_topaz         No 
do_topaz_filaments         No 
do_topaz_pick         No 
do_topaz_train         No 
do_topaz_train_parts         No 
do_write_fom_maps         No 
fn_input_autopick Schemes/amyprep/ctffind/micrographs_ctf.star 
fn_ref3d_autopick         "" 
fn_refs_autopick         "" 
fn_topaz_exe relion_python_topaz 
   gpu_ids         "" 
helical_nr_asu          3 
helical_rise       4.75 
helical_tube_kappa_max       0.07 
helical_tube_length_min        500 
helical_tube_outer_diameter        200 
  highpass         -1 
log_adjust_thr          0 
log_diam_max        250 
log_diam_min        200 
log_invert         No 
log_maxres         20 
log_upper_thr        999 
   lowpass         20 
maxstddevnoise_autopick        1.1 
min_dedicated        112 
minavgnoise_autopick       -999 
mindist_autopick        100 
    nr_mpi          6 
nr_threads         18 
other_args         "" 
psi_sampling_autopick          5 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-cpu-devel.sh 
 queuename    openmpi 
ref3d_sampling "30 degrees" 
ref3d_symmetry         C1 
    shrink          0 
threshold_autopick       0.05 
topaz_filament_threshold         -5 
topaz_hough_length         -1 
topaz_model         "" 
topaz_nr_particles         -1 
topaz_other_args         "" 
topaz_particle_diameter         -1 
topaz_train_parts         "" 
topaz_train_picks         "" 
   use_gpu         No 
 
