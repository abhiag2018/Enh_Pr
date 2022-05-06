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

PYTHONPATH=NN-feature-prep/bin:$PYTHONPATH
# base_dir="/projects/li-lab/agarwa/CUBE/DeepTact/DeepTact-Pipeline/storeDir"
# code_dir="/projects/li-lab/agarwa/CUBE/DeepTact/DeepTact-Pipeline/DeepTact_NF_code/bin"
base_dir="."
code_dir="."

job=$1; shift
cell=$1; shift
num_rep=$1; shift
eval_cell=$1; shift

# if [ $job == "train" ] && [ "$#" -ne 2 ]; then
#     echo "2 arguments needed : train <celltype> <numRep>"
# fi

# if [ $job == "split" ] && [ "$#" -ne 2 ]; then
#     echo "split: 2 arguments needed : split <celltype> <numRep>"
# fi

# if [ $job == "test" ] && [ "$#" -lt 5 ]; then
#     echo ">=5 arguments needed : test <celltype> <numRep> <evalCell> <appendStr> <test-file> <bootstrap1> <bootstrap2>.."
# fi

#python -c "import pandas as pd; df=pd.read_csv('${base_dir}/${cell}/P-E/_train.csv'); df.data=df.data.apply(lambda x:'data-train/'+x); df.to_csv('${base_dir}/${cell}/P-E/train.csv',index=False)"
#python -c "import pandas as pd; df=pd.read_csv('${base_dir}/${cell}/P-E/_val.csv'); df.data=df.data.apply(lambda x:'data-val/'+x); df.to_csv('${base_dir}/${cell}/P-E/val.csv',index=False)"

python $code_dir/DeepTact_5.py $job ${base_dir}/${cell} ${num_rep} ${eval_cell} "$@"
#python $code_dir/DeepTact_2.py $job ${base_dir}/${cell} P-E ${num_rep} ${eval_cell} "$@"

#sbatch scriptGPU.sh train mouse_cube/rep1/type77 3
#sbatch scriptGPU.sh train mouse_cube_type_test 3
