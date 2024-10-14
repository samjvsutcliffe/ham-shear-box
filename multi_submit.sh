#!/bin/bash
module load gcc
module load mvapich2
module load aocl
export MV2_ENABLE_AFFINITY=0
#sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit
rm -r output-*


for ref in 4
do
    #for l in 100000 125000 150000 175000 200000 225000 250000 275000 300000
    for l in 100000 200000 300000
    do
        export REFINE=$ref
        export LOAD=$l
        sbatch batch_shear-box.sh
    done
done


#export REFINE=1
#sbatch batch_shear-box.sh 
#export REFINE=2
#sbatch batch_shear-box.sh 
#export REFINE=4
#sbatch batch_shear-box.sh 
#export REFINE=8
#sbatch batch_shear-box.sh 
#export REFINE=16
#sbatch batch_shear-box.sh 
#export REFINE=32
#sbatch batch_shear-box.sh 
#export REFINE=64
#sbatch batch_shear-box.sh 
