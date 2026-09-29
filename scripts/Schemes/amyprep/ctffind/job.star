
# version 50001

data_job

_rlnJobTypeLabel             relion.ctffind.ctffind4
_rlnJobIsContinue                       0
_rlnJobIsTomo                           0


# version 50001

data_joboptions_values

loop_ 
_rlnJobOptionVariable #1 
_rlnJobOptionValue #2 
       box        512 
   ctf_win         -1 
      dast        100 
     dfmax      50000 
     dfmin       5000 
    dfstep        500 
do_phaseshift         No 
  do_queue        No
fn_ctffind_exe /public/EM/ctffind/ctffind.exe 
input_star_mics Schemes/amyprep/motioncorr/corrected_micrographs.star
min_dedicated     112 
    nr_mpi         24 
other_args         "" 
 phase_max        180 
 phase_min          0 
phase_step         10 
      qsub     sbatch 
qsubscript /public/EM/RELION/relion-slurm-cpu-devel.sh
 queuename    openmpi 
    resmax          5 
    resmin         30 
slow_search         No 
use_given_ps        Yes 
  use_noDW         No 
