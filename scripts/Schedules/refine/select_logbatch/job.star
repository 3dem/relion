
# version 30001

data_job

_rlnJobTypeLabel             Select
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
discard_label rlnImageName 
discard_sigma          4 
do_class_ranker        Yes 
do_discard         No 
  do_queue         No 
 do_random         No 
do_recenter         No 
do_regroup         No 
do_remove_duplicates         No 
do_select_values         No 
  do_split         No 
duplicate_threshold         30 
 fn_coords         "" 
   fn_data         "" 
    fn_mic         "" 
  fn_model Schedules/refine/class2d_logbatch/run_it025_optimiser.star 
image_angpix         -1 
min_dedicated         24 
 nr_groups          1 
  nr_split         -1 
other_args         "" 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
rank_threshold       0.35 
select_label rlnCtfMaxResolution 
select_maxval      9999. 
select_minval     -9999. 
split_size        100 
 