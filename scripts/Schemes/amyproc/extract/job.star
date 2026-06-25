
# version 50001

data_job

_rlnJobTypeLabel             relion.extract.helical
_rlnJobIsContinue                       0
_rlnJobIsTomo                           0
 

# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
bg_diameter         -1 
black_dust         -1 
coords_suffix Schemes/amyproc/autopick_trace/autopick.star 
do_extract_helix        Yes 
do_float16        Yes 
do_fom_threshold         No 
 do_invert        Yes 
   do_norm        Yes 
  do_queue         No 
do_recenter         No 
do_reextract         No 
do_rescale        Yes 
do_reset_offsets         No 
extract_size        768 
fndata_reextract ""
helical_bimodal_angular_priors        Yes 
helical_nr_asu          3 
helical_rise       4.75 
helical_tube_outer_diameter        200 
min_dedicated         24 
minimum_pick_fom          0 
    nr_mpi         12 
other_args         "" 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-cpu-devel.sh 
 queuename    openmpi 
recenter_x          0 
recenter_y          0 
recenter_z          0 
   rescale        128 
 star_mics $$input_mics 
white_dust         -1 
 
