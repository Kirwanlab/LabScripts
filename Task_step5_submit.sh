#!/bin/bash





workDir=~/compute/STT_reml					###??? update this
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/TS5_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR5.txt \
-e ${outDir}/error_sttR5.txt \
Task_step5_meanBetas.sh
