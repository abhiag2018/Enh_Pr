#!/bin/bash
#SBATCH --mem=1G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###

basedir="../nxf-scripts"

inp_config1="nextflow.config"
inp_config2="conf/pr-enh-prep.config"
refgen="mm10"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	run $basedir/pr-enh-prep.nf \
	--refgen $refgen \
	-profile slurm \
	-w "./work" -with-timeline \
	"$@"
	# -resume
# --dev \
