#!/bin/bash


# Written by Nathan Muncy on 2/1/18

# Wrapper script for the sbatch job


### Notes
#
# 1) This will create output txt files for each sbatch submission in $outDir
#
# 2) output files will be organized by type and time of job submission



# Set up output dirs, good for troubleshooting
workDir=~/compute/<something>

slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/ashs_${time}

mkdir -p $outDir


# Do work
cd $workDir

for i in s*; do

    sbatch -o ${outDir}/output_ashstest.txt \
    -e ${outDir}/error_ashstest.txt \
    ashs_step1_sbatch.sh $workDir $i

    sleep 1
done
