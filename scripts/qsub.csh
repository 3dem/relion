#!/bin/tcsh
#$ -pe XXXqueueXXX XXXmpinodesXXX
#$ -l dedicated=XXXthreadsXXX 
#$ -e XXXerrfileXXX
#$ -o XXXoutfileXXX
#$ -cwd
#$ -S /bin/tcsh

# Environment
source ~/.cshrc

mpiexec -n XXXmpinodesXXX  XXXcommandXXX
