#!/bin/bash
#SBATCH -q training 
#SBATCH -p gpu
#SBATCH --time=11-01:00:00
#SBATCH --gres=gpu:v100:1
#SBATCH --mem=32G

##SBATCH -q dev 
##SBATCH --time=8:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate /projects/li-lab/agarwa/conda_envs/NN-gpu

base_dir="."
code_dir="../nxf_TestSet/nxf-scripts/bin/"
PYTHONPATH=$code_dir:$PYTHONPATH

job=$1; shift
cell=$1; shift
num_rep=$1; shift
eval_cell=$1; shift
model_id=$1; shift

python $code_dir/DeepTact_6.py $job ${base_dir}/${cell} ${num_rep} ${eval_cell} ${model_id} "$@"

#sbatch scriptGPU.sh train type74 3
#sbatch scriptGPU.sh test type74/data/dataset-test.h5 type74 3 type74_test 20220512-134252-XPD 
