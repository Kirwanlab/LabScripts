#!/bin/bash




### --- Notes
#
# 1) this script will construct T1, T2, EPI data and organize output
#		according to BIDS formatting
#
# 2) written so you can just update $subjList and rerun the whole script



###??? change these variables/arrays
rawDir=/Volumes/Yorick/MriRawData							# location of raw data
workDir=/Volumes/Yorick/Nate_work/STT_reml					# desired working directory

subjList=(1295 1859 1875 1949)								# list of subjects
session=STT													# scanning session - for raw data organization (ses-STT)
task=task-stt												# name of task, for epi data naming

epiDirs=(STUDY TEST1{1,2} TEST2{1..4})						# epi dicom directory name/prefix
t1Dir=t1													# t1 ditto
t2Dir=HHR													# t2 ditto
blipDir=()													# blip ditto - one per scanning phase



### set up BIDS parent dirs
for i in derivatives sourcedata stimuli; do
	if [ ! -d ${workDir}/$i ]; then
		mkdir -p ${workDir}/$i
	fi
done



for i in ${subjList[@]}; do

	### set up BIDS data dirs
	anatDir=${workDir}/rawdata/sub-${i}/anat
	funcDir=${workDir}/rawdata/sub-${i}/func
	derivDir=${workDir}/derivatives/sub-${i}

	if [ ! -d $anatDir ]; then
		mkdir -p $anatDir $funcDir $derivDir
	fi


	### construct data
	dataDir=${rawDir}/sub-${i}/ses-${session}/dicom

	# t1 data
	if [ ! -f ${anatDir}/sub-${i}_T1w.nii.gz ]; then
		dcm2niix -b y -ba y -z y -o $anatDir -f sub-${i}_T1w ${dataDir}/${t1Dir}*/
	fi


	# t2
	if [ ! -f ${anatDir}/sub-${i}_T2w.nii.gz ]; then
		dcm2niix -b y -ba y -z y -o $anatDir -f sub-${i}_T2w ${dataDir}/${t2Dir}*/
	fi


	# epi
	for j in ${!epiDirs[@]}; do
		pos=$(($j+1))
		if [ ! -f ${funcDir}/sub-${i}_${task}_run-${pos}_bold.nii.gz ]; then
			dcm2niix -b y -ba y -z y -o $funcDir -f sub-${i}_${task}_run-${pos}_bold ${dataDir}/${epiDirs[$j]}*/
		fi
	done


	# blip
	if [ ! -z $blipDir ]; then
		for k in ${!blipDir[@]}; do
			pos=$(($k+1))
			if [ ! -f ${funcDir}/sub-${i}_phase-${pos}_blip.nii.gz ]; then
				dcm2niix -b y -ba y -z y -o $funcDir -f sub-${i}_phase-${pos}_blip ${dataDir}/${blipDir}*/
			fi
		done
	fi
done
