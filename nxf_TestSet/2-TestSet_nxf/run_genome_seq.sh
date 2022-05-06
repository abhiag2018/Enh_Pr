#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###

# configdir="$basedir/config_hg38_mm10"
# cp $configdir/genome-seq-prep.config .
# cp $configdir/pr-enh-prep.config .
# cp $basedir/nextflow-1.config nextflow.config

refgen="mm10"
inp_config1="nextflow.config"
inp_config2="conf/genome-seq-prep.config"
inp_config3="conf/pr-enh-prep.config"

basedir="../nxf-scripts"
nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/genome-seq-prep.nf \
	--refgen $refgen \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	"$@"
	# -resume
	# --dev 
