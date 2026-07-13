#!/bin/bash
module load aocc/5.0.0
module load aocl/5.0.0
module load mvapich2
export MV2_ENABLE_AFFINITY=0
sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit
#rm -r output-*

# 8 16
for ref in 8
do
    for l in 100000 200000 300000
    #for l in 200000
    #for l in 25000 50000 75000 100000 125000 150000 175000 200000 225000 250000 275000 300000
    do
        export REFINE=$ref
        export LOAD=$l
        sbatch sb-d-mc.sh
    done
done

#for ref in 8 16
#do
#    #for l in 50000 100000 125000 150000 175000 200000 225000 250000 275000 300000
#    #for l in 125000 150000 175000 225000 250000 275000
#    #for l in 100000 150000 200000 250000 300000 350000
#    for l in 100000 200000 300000
#    do
#        export REFINE=$ref
#        export LOAD=$l
#        sbatch sb-d-mc-big.sh
#    done
#done
#
#for ref in 32
#do
#    #for l in 50000 100000 125000 150000 175000 200000 225000 250000 275000 300000
#    #for l in 125000 150000 175000 225000 250000 275000
#    #for l in 100000 150000 200000 250000 300000 350000
#    for l in 100000 200000 300000
#    do
#        export REFINE=$ref
#        export LOAD=$l
#        sbatch sb-d-mc-verybig.sh
#    done
#done


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
##sbatch batch_shear-box.sh 
