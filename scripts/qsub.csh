#!/bin/tcsh
#$ -pe XXXqueueXXX XXXmpinodesXXX
#$ -l dedicated=XXXthreadsXXX 
#$ -e XXXerrfileXXX
#$ -o XXXoutfileXXX
#$ -cwd
#$ -S /bin/tcsh

# Environment
source ~/.cshrc

mpiexec -mca orte_forward_job_control 1 -n XXXmpinodesXXX  XXXcommandXXX
