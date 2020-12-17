
# version 30001

data_job

_rlnJobTypeLabel             JoinStar
_rlnJobIsContinue                       0
 

# version 30001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
    do_mic         No 
    do_mov         No 
   do_part        Yes 
  do_queue         No 
   fn_mic1         "" 
   fn_mic2         "" 
   fn_mic3         "" 
   fn_mic4         "" 
   fn_mov1         "" 
   fn_mov2         "" 
   fn_mov3         "" 
   fn_mov4         "" 
  fn_part1 $$all_particles 
  fn_part2 Schedules/class2d/class2d_rest/run_it250_data.star 
  fn_part3         "" 
  fn_part4         "" 
min_dedicated         24 
other_args         "" 
      qsub       qsub 
qsubscript /public/EM/RELION/relion/bin/relion_qsub.csh 
 queuename    openmpi 
 
