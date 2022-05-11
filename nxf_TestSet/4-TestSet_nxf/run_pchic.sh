#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00


###
NXF_OPTS='-Xms16g -Xmx64g'

basedir="../nxf-scripts"

refgen="mm10"
inp_config1="nextflow.config"
inp_config2="conf/pchic-prep.config"
inp_config3="conf/pr-enh-prep.config"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/pchic-prep.nf \
	-profile slurm \
	-w "./work" -with-timeline \
	--refgen $refgen \
	"$@"

#old # ./run_pchic.sh --mouse_sample ORSAM17820-2 
#new ./run_pchic.sh
## 	--dev
