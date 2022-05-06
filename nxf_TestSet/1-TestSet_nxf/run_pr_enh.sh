#!/bin/bash
#SBATCH --mem=1G
#SBATCH --qos batch 
#SBATCH --time 6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###

basedir="../nxf-scripts"
# configdir="$basedir/config_hg38_mm10"
# cp $configdir/pr-enh-prep.config .
# cp $basedir/nextflow-1.config nextflow.config


inp_config1="nextflow.config"
inp_config2="conf/pr-enh-prep.config"
refgen="mm10"

basedir="/projects/li-lab/agarwa/CUBE/DeepTact/DeepTact-Pipeline/DeepTact_NF_code"
nextflow -C $inp_config1 \
	-C $inp_config2 \
	run $basedir/pr-enh-prep.nf \
	--refgen $refgen \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	"$@"
	# -resume
# --dev \
