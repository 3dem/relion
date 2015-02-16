#!/bin/bash

BINDIR=$(cd ../.. ; echo $PWD)/bin
echo "...using binaries in $BINDIR"
#ls $BINDIR

shopt -s expand_aliases
alias refine="$BINDIR/relion_refine"

refine \
 --o Class3D/run1\
 --i particles_autopick_sort_class2d.star\
 --particle_diameter 200\
 --angpix 1\
 --ref 3i3e_lp50A.mrc\
 --firstiter_cc\
 --ini_high 50\
 --ctf\
 --iter 25\
 --tau2_fudge 2\
 --K 4\
 --flatten_solvent\
 --zero_mask\
 --oversampling 1\
 --healpix_order 2\
 --offset_range 5\
 --offset_step 2\
 --sym C1\
 --norm\
 --scale\
 --j 1\
 --memory_per_thread 4


