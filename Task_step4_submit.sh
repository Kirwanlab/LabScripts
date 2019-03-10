#!/bin/bash





workDir=~/compute/STT_reml
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR4_${time}

mkdir -p $outDir


sbatch \
-o ${outDir}/output_sttR4.txt \
-e ${outDir}/error_sttR4.txt \
Task_step4_grpAnalysis.sh
