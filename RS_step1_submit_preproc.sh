#!/bin/bash


# written by Nathan Muncy on 12/3/19


parDir=~/compute/MNI_NREM/MRI_processed			###??? update this
rawDir=${parDir}/rawdata
workDir=${parDir}/derivatives
scriptDir=${parDir}/code

phaseArr=(Sleep Wake)							###??? phases of experiment, corresponds to sessions (ses-Sleep, ses-Wake). Required

slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/RS1_${time}

mkdir -p $outDir


cd $rawDir

for i in s*; do
	for j in ${phaseArr[@]}; do
		if [ ! -f ${workDir}/${i}/ses-${j}/Resting_${j}_scale+tlrc.HEAD ]; then

		    sbatch \
		    -o ${outDir}/output_RS1_${i}_${j}.txt \
		    -e ${outDir}/error_RS1_${i}_${j}.txt \
		    ${scriptDir}/RS_step1_sbatch_preproc.sh $i $j

		    sleep 1
		fi
	done
done
