#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 2-12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

basedir="../nxf-scripts"
refgen="mm10"

inp_config1="nextflow.config"
inp_config2="conf/feature-prep.config"
inp_config3="conf/pr-enh-prep.config"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/NN-feature-prep.nf \
	-profile slurm \
	-w "./work" -with-timeline \
	--refgen $refgen \
	"$@"
	# -resume ${resumeID}
	## --dev

# ./run_feature_gen.sh -resume
