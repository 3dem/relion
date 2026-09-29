
# version 50001

data_job

_rlnJobTypeLabel             relion.import.movies
_rlnJobIsContinue                       0
_rlnJobIsTomo                           0
 

# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
        Cs        2.7 
        Q0        0.1 
    angpix      0.824 
beamtilt_x          0 
beamtilt_y          0 
  do_other         No 
  do_queue         No 
    do_raw        Yes 
fn_in_other    ref.mrc 
 fn_in_raw Tiff/*/*/*.tiff 
    fn_mtf         "" 
is_multiframe        Yes 
        kV      300.0 
min_dedicated         24 
 node_type "Particle coordinates (*.box, *_pick.star)" 
optics_group_name opticsGroup1 
optics_group_particles         "" 
other_args " --do_at_most $$do_at_most" 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-cpu-devel.sh 
 queuename    openmpi 
 
