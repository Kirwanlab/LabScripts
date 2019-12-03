#!/bin/bash


# written by Nathan Muncy on 12/3/19


parDir=~/compute/MNI_NREM/MRI_processed			###??? update this
rawDir=${parDir}/rawdata
workDir=${parDir}/derivatives
scriptDir=${parDir}/code

phaseArr=(Sleep Wake)							###??? update this

slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/RS2_${time}

mkdir -p $outDir


cd $rawDir

for i in s*; do
	for phase in ${phaseArr[@]}; do
		if [ ! -f ${workDir}/${i}/ses-${phase}/errts_fanaticor+tlrc.HEAD ]; then

		    sbatch \
		    -o ${outDir}/output_RS2_${i}_${phase}.txt \
		    -e ${outDir}/error_RS2_${i}_${phase}.txt \
		    ${scriptDir}/RS_step2_regress.sh $i $phase

		    sleep 1
		fi
	done
done
