#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 2-12:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'

alias hic="cd ../pchic_mouse_cube/"
alias wd="cd `pwd`"
basedir="../nxf-scripts"
refgen="mm10"
#chunks=11

inp_config1="nextflow.config"
inp_config2="conf/feature-prep.config"
inp_config3="conf/pr-enh-prep.config"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/NN-feature-prep.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--refgen $refgen \
	"$@"
	# -resume ${resumeID}
	## --dev

# ./run_feature_gen.sh -resume
