#!/bin/bash


# written by Nathan Muncy on 11/7/17


### Wrapper script for RestingState_step1_sbatch_preproc.sh


workDir=~/compute/MNI_RDoC/MRI_processed
scriptDir=${workDir}/Scripts

# Create output dir for trouble shooting
slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/resting_${time}
mkdir -p $outDir


cd $workDir

for i in P*; do

    sbatch \
    -o ${outDir}/output_resting_${i}.txt \
    -e ${outDir}/error_resting_${i}.txt \
    ${scriptDir}/RestingState_step1_sbatch_preproc.sh $workDir $i

    sleep 1
done
