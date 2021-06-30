
# version 30001

data_job

_rlnJobTypeLabel             relion.autopick
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
    angpix         -1 
angpix_ref         -1 
do_amyloid         No 
do_ctf_autopick        Yes 
do_ignore_first_ctfpeak_autopick         No 
do_invert_refs        Yes 
    do_log         $$do_log
do_pick_helical_segments         No 
  do_queue         No 
do_read_fom_maps         No 
  do_ref3d        Yes 
   do_refs         No 
  do_topaz        $$do_topaz
do_topaz_pick        Yes 
do_topaz_train         No 
do_topaz_train_parts         No 
do_write_fom_maps         No 
fn_input_autopick Schemes/proc/select_mics/micrographs.star
fn_ref3d_autopick         "" 
fn_refs_autopick         "" 
   gpu_ids          0
helical_nr_asu          1 
helical_rise         -1 
helical_tube_kappa_max        0.1 
helical_tube_length_min         -1 
helical_tube_outer_diameter        200 
  highpass         -1 
log_adjust_thr          0 
log_diam_max        180 
log_diam_min        150 
log_invert         No 
log_maxres         20 
log_upper_thr        999 
   lowpass         20 
maxstddevnoise_autopick         "" 
min_dedicated         24 
minavgnoise_autopick         "" 
mindist_autopick         "" 
    nr_mpi          2 
other_args         "" 
psi_sampling_autopick          5 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
ref3d_sampling "30 degrees" 
ref3d_symmetry         "" 
    shrink          0 
threshold_autopick         "" 
topaz_model $$topaz_model 
topaz_nr_particles        900 
topaz_other_args         "" 
topaz_particle_diameter      100.0 
topaz_train_parts         "" 
topaz_train_picks         "" 
   use_gpu        Yes 
 
