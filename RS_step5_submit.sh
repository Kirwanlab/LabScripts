#!/bin/bash


# written by Nathan Muncy on 12/3/19


parDir=~/compute/MNI_NREM/MRI_processed			###??? update this
workDir=${parDir}/derivatives
scriptDir=${parDir}/code
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/RS5_${time}

mkdir -p $outDir

sbatch \
-o ${outDir}/output_RS5.txt \
-e ${outDir}/error_RS5.txt \
${scriptDir}/RS_step5_grpAnalysis.sh

