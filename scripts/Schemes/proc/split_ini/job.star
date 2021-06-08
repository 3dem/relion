
# version 30001

data_job

_rlnJobTypeLabel             relion.select
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
discard_label rlnImageName 
discard_sigma          4 
do_class_ranker         No 
do_discard         No 
  do_queue         No 
 do_random         No 
do_recenter         No 
do_regroup         No 
do_remove_duplicates         No 
do_select_values         No 
  do_split        Yes 
duplicate_threshold         30 
 fn_coords         "" 
   fn_data Schemes/proc/extract_ini/particles.star 
    fn_mic         "" 
  fn_model         "" 
image_angpix         -1 
min_dedicated         24 
 nr_groups          1 
  nr_split         -1 
other_args         "" 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
rank_threshold        0.5 
select_label rlnCtfFigureOfMerit 
select_maxval      9999. 
select_minval     -9999. 
split_size      10000 
 
