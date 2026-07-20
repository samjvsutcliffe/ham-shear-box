#!/bin/bash
module load aocc/5.0.0
module load aocl/5.0.0
module load mvapich2
export MV2_ENABLE_AFFINITY=0
sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit
#rm -r output-*

# 0.999
for d in 0 0.5 0.9 0.99 0.999
do
    for o in 1
    do
        for ref in 4
        do
            #for l in 100000 125000 150000 175000 200000 225000 250000 275000 300000
            #for l in 125000 150000 175000 225000 250000 275000
            for l in 100000 200000 300000
            do
                export OVER=$o
                export REFINE=$ref
                export LOAD=$l
                export DAMAGE=$d
                sbatch sb-pdr.sh
            done
        done
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
