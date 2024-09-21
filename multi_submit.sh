module load gcc
module load mvapich2
module load aocl
export MV2_ENABLE_AFFINITY=0
#sbcl --dynamic-space-size 16000  --disable-debugger --load "build_step.lisp" --quit
rm -r output-*

export REFINE=4
export LOAD=30000
sbatch batch_shear-box.sh 
export LOAD=70000
sbatch batch_shear-box.sh 
export LOAD=100000
sbatch batch_shear-box.sh 
export LOAD=150000
sbatch batch_shear-box.sh 
export LOAD=200000
sbatch batch_shear-box.sh 
#export LOAD=300000
#sbatch batch_shear-box.sh 
#
#export LOAD=100000

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
