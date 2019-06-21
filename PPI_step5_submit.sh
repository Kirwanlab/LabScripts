#!/bin/bash





workDir=~/compute/AutismOlfactory
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/PPI5_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_PPI5.txt \
-e ${outDir}/error_PPI5.txt \
PPI_step5_grpAnalysis.sh
