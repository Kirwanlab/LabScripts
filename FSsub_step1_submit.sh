#!/bin/bash


# Written by Nathan Muncy on 2/6/18


workDir=~/compute/<something>
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/FSsub_${time}

mkdir -p $outDir


cd $workDir

for i in s*; do

    sbatch \
    -o ${outDir}/output_FSsub_${i}.txt \
    -e ${outDir}/error_FSsub_${i}.txt \
    ${scriptDir}/FSsub_step1_segment.sh $i $workDir

    sleep 1
done


