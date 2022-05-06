#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 2-6:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate cube

###
NXF_OPTS='-Xms16g -Xmx64g'


basedir="../nxf-scripts"
configdir="$basedir/config_hg38_mm10"
# cp $configdir/co-score-prep.config .
# cp $configdir/pr-enh-prep.config .
# cp $basedir/nextflow-1.config nextflow.config

refgen="mm10"
inp_config1="nextflow.config"
inp_config2="conf/co-score-prep.config"
inp_config3="conf/pr-enh-prep.config"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/co-score-prep.nf \
	-profile slurm -w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--refgen $refgen \
	"$@"
	# -resume \
	# --dev
