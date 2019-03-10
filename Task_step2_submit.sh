#!/bin/bash




###??? update these
workDir=~/compute/STT_reml
scriptDir=${workDir}/code
slurmDir=${workDir}/derivatives/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/sttR2_${time}

mkdir -p $outDir

cd ${workDir}/derivatives
for i in sub*; do

	[ $i == sub-1295 ]; test=$?

    sbatch \
    -o ${outDir}/output_sttN2_${i}.txt \
    -e ${outDir}/error_sttN2_${i}.txt \
    ${scriptDir}/Task_step2_sbatch_regress.sh $i $test

    sleep 1
done
