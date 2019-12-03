#!/bin/bash


# written by Nathan Muncy on 12/3/19


parDir=~/compute/MNI_NREM/MRI_processed			###??? update this
workDir=${parDir}/derivatives
scriptDir=${parDir}/code

phaseArr=(Sleep Wake)							###??? update this

slurmDir=${workDir}/Slurm_out
time=`date '+%Y_%m_%d-%H_%M_%S'`
outDir=${slurmDir}/RS3_${time}

mkdir -p $outDir


cd ${parDir}/rawdata

for i in s*; do
	for j in ${phaseArr[@]}; do
		if [ ! -f ${workDir}/${i}/ses-${j}/FINAL_Seed_lLing+tlrc.HEAD ]; then		###??? Name of a seed here from job script

		    sbatch \
		    -o ${outDir}/output_RS3_${i}_${j}.txt \
		    -e ${outDir}/error_RS3_${i}_${j}.txt \
		    ${scriptDir}/RS_step3_seedTS.sh $i $j

		    sleep 1
		fi
	done
done
