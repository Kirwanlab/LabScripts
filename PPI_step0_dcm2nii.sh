#!/bin/bash




### --- Notes
#
# 1) this script will construct T1, T2, EPI data and organize output
#		according to BIDS formatting
#
# 2) written so you can just update $subjList and rerun the whole script



###??? change these variables/arrays
rawDir=/Volumes/Yorick/Nate_work/De_Identified_BO_Data							# location of raw data
workDir=/Volumes/Yorick/Nate_work/AutismOlfactory					# desired working directory

task=task-AO												# name of task, for epi data naming
epiDirs=(BO_Run{1..3})										# epi dicom directory name/prefix
t1Dir=t1													# t1 ditto




### set up BIDS parent dirs
for i in derivatives sourcedata stimuli; do
	if [ ! -d ${workDir}/$i ]; then
		mkdir -p ${workDir}/$i
	fi
done


cd $rawDir
for i in S*; do


	### set up BIDS data dirs
	subj=${i##*_}
	anatDir=${workDir}/rawdata/sub-${subj}/anat
	funcDir=${workDir}/rawdata/sub-${subj}/func
	derivDir=${workDir}/derivatives/sub-${subj}

	if [ ! -d $anatDir ]; then
		mkdir -p $anatDir $funcDir $derivDir
	fi


	### construct data
	dataDir=${rawDir}/${i}/Res*

	# t1 data
	if [ ! -f ${anatDir}/sub-${subj}_T1w.nii.gz ]; then
		dcm2niix -b y -ba y -z y -o $anatDir -f sub-${i}_T1w ${dataDir}/${t1Dir}*/
	fi


	# epi
	for j in ${!epiDirs[@]}; do
		pos=$(($j+1))
		if [ ! -f ${funcDir}/sub-${subj}_${task}_run-${pos}_bold.nii.gz ]; then
			dcm2niix -b y -ba y -z y -o $funcDir -f sub-${subj}_${task}_run-${pos}_bold ${dataDir}/${epiDirs[$j]}*/
		fi
	done
done
