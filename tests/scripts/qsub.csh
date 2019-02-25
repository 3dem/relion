#!/bin/tcsh
#$ -pe XXXqueueXXX XXXnodesXXX
#$ -l dedicated=XXXdedicatedXXX 
#$ -e XXXerrfileXXX
#$ -o XXXoutfileXXX
#$ -A Relion 
#$ -cwd
#$ -S /bin/tcsh

mpiexec -mca orte_forward_job_control 1 -n XXXmpinodesXXX  XXXcommandXXX
