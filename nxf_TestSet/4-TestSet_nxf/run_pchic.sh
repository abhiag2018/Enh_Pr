#!/bin/bash
#SBATCH --mem=10G
#SBATCH --qos batch 
#SBATCH --time 3-00:00:00


###
NXF_OPTS='-Xms16g -Xmx64g'

# basedir="/projects/li-lab/agarwa/CUBE/DeepTact/DeepTact-Pipeline/DeepTact_NF_code"
basedir="../nxf-scripts"
configdir="$basedir/config_mm10"
# cp $configdir/pchic-prep_mm.config .
# cp $basedir/nextflow-1.config nextflow.config

refgen="mm10"
inp_config1="nextflow.config"
inp_config2="conf/pchic-prep.config"
inp_config3="conf/pr-enh-prep.config"

nextflow -C $inp_config1 \
	-C $inp_config2 \
	-C $inp_config3 \
	run $basedir/pchic-prep.nf \
	-profile slurm \
	-w "/fastscratch/agarwa/nf-tmp/work" -with-timeline \
	--refgen $refgen \
	"$@"

#old # ./run_pchic.sh --mouse_sample ORSAM17820-2 
#new ./run_pchic.sh
## 	--dev
