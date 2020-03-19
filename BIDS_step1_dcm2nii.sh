#!/bin/bash





### --- Notes
#
# 1) this script will construct T1, T2, EPI data and organize output
#		according to BIDS formatting
#
# 2) written so you can just update $subjList and rerun the whole script
#
# 3) Will require input for constructing dataset_description.json
#
# 4) Make sure to update the jq task input where the epi Json files are being input!
#
# 5) If an experiment has two phases (Study, Test), where Study has 2 runs and Test 4,
#		then: expPhaseRuns=(2 4)
#		and Task runs 1-2 are Study, 3-6 are Test.




###??? change these variables/arrays
rawDir=/Volumes/Yorick/MriRawData							# location of raw data
workDir=/Volumes/Yorick/Temporal							# desired working directory
tempDir=/Volumes/Yorick/Templates

subjList=(`cat ${workDir}/Experiment3/subjList.txt`)		# a column of subjects (e.g. sub-111 sub-112 sub-113)
session=Temporal											# scanning session - for raw data organization (ses-STT)
task=TemporalMST											# name of task, for epi data naming

epiDirs=(Temporal_B{1..3})									# epi dicom directory name/prefix
t1Dir=t1													# t1 ditto
t2Dir=HHR													# t2 ditto

blipDir=(Reverse_Blip)										# blip ditto - names one per scanning phase
blipEnc=PA 													# alphanumeric string for blip direction
expPhaseRuns=(3)											# Note 5




### Check for jq
# jq will be used to append Json files

which jq >/dev/null 2>&1

if [ $? != 0 ]; then
	echo >&2
	echo "Software jq is required: download from https://stedolan.github.io/jq/ and add it to your \$PATH. Exit 1" >&2
	echo >&2
	exit 1
fi




### Set up BIDS parent dirs
for i in derivatives sourcedata stimuli rawdata; do
	if [ ! -d ${workDir}/$i ]; then
		mkdir -p ${workDir}/$i
	fi
done




### Write dataset_description.json
# This will request input

if [ ! -s ${workDir}/rawdata/dataset_description.json ]; then

	echo -e "\nNote: title below must be supplied in quotations"
	echo -e "\te.g. \"This is my title\"\n"
	read -p 'Please enter title of the manuscript:' title

	echo -e "\n\nNote: authors must be within quotes and separated by a comma & space"
	echo -e "\te.g. \"Nate M. Muncy\", \"C. Brock Kirwan\"\n"
	read -p 'Please enter authors:' authors

	cat > ${workDir}/rawdata/dataset_description.json << EOF
{
	"Name": $title,
	"BIDSVersion": "1.1.1",
	"License": "CCo",
	"Authors": [$authors]
}
EOF
fi



### Work
for i in ${subjList[@]}; do

	### set up BIDS data dirs
	anatDir=${workDir}/rawdata/${i}/anat
	funcDir=${workDir}/rawdata/${i}/func
	fmapDir=${workDir}/rawdata/${i}/fmap
	derivDir=${workDir}/derivatives/${i}

	for j in {anat,func,fmap,deriv}Dir; do
		if [ ! -d $(eval echo \${$j}) ]; then
			mkdir -p $(eval echo \${$j})
		fi
	done


	### construct data
	dataDir=${rawDir}/${i}/ses-${session}/dicom

	# t1 data
	if [ ! -f ${anatDir}/${i}_T1w.nii.gz ]; then

		dcm2niix -b y -ba y -z y -o $anatDir -f tmp_${i}_T1w ${dataDir}/${t1Dir}*/

		# Deface by Brock Kirwan
		3dAllineate -base ${anatDir}/tmp_${i}_T1w.nii.gz -input ${tempDir}/mean_reg2mean.nii.gz -prefix ${anatDir}/tmp_mean_reg2mean_aligned.nii -1Dmatrix_save ${anatDir}/tmp_allineate_matrix
		3dAllineate -base ${anatDir}/tmp_${i}_T1w.nii.gz -input ${tempDir}/facemask.nii.gz -prefix ${anatDir}/tmp_facemask_aligned.nii -1Dmatrix_apply ${anatDir}/tmp_allineate_matrix.aff12.1D
		3dcalc -a ${anatDir}/tmp_facemask_aligned.nii -b ${anatDir}/tmp_${i}_T1w.nii.gz -prefix ${anatDir}/${i}_T1w.nii.gz -expr "step(a)*b"
		mv ${anatDir}/tmp_${i}_T1w.json ${anatDir}/${i}_T1w.json
		rm ${anatDir}/tmp*
	fi


	# t2
	if [ ! -f ${anatDir}/${i}_T2w.nii.gz ]; then
		dcm2niix -b y -ba y -z y -o $anatDir -f ${i}_T2w ${dataDir}/${t2Dir}*/
	fi


	# epi
	for j in ${!epiDirs[@]}; do

		pos=$(($j+1))

		if [ ! -f ${funcDir}/${i}_task-${task}_run-${pos}_bold.nii.gz ]; then
			dcm2niix -b y -ba y -z y -o $funcDir -f ${i}_task-${task}_run-${pos}_bold ${dataDir}/${epiDirs[$j]}*/
		fi

		# Json append by Brock Kirwan
		# updated by Nate to receive argument
		funcJson=${funcDir}/${i}_task-${task}_run-${pos}_bold.json
		taskExist=$(cat $funcJson | jq '.TaskName')
		if [ "$taskExist" == "null" ]; then
			jq --arg holdTask "$task" '. |= . + {"TaskName":$holdTask}' $funcJson > ${derivDir}/tasknameadd.json
			rm $funcJson && mv ${derivDir}/tasknameadd.json $funcJson
		fi
	done


	# blip
	if [ ! -z $blipDir ]; then

		runningC=0
		for k in ${!blipDir[@]}; do

			pos=$(($k+1))
			if [ ! -f ${fmapDir}/${i}_ses-${session}_dir-${blipEnc}_run-${pos}_phasediff.nii.gz ]; then
				dcm2niix -b y -ba y -z y -o $fmapDir -f ${i}_ses-${session}_dir-${blipEnc}_run-${pos}_phasediff ${dataDir}/${blipDir}*/
			fi

			## append Json for intended use
			fmapJson=${fmapDir}/${i}_ses-${session}_dir-${blipEnc}_run-${pos}_phasediff.json
			intendExist=$(cat $fmapJson | jq '.IntendedFor')
			if [ "$intendExist" == "null" ]; then

				# build list of phase runs
				runningC=$(($runningC+${expPhaseRuns[$k]}))
				unset subjArr
				run=1; while [ $run -le $runningC ]; do
					arrPos=$(($run-1))
					subjArr[$arrPos]=func/${i}_task-${task}_run-${run}_bold.nii.gz
					let run+=1
				done

				printf '%s\n' "${subjArr[@]}" | jq -R . > ${derivDir}/tmpSubj.json
				jq --slurpfile jsonSubjArr ${derivDir}/tmpSubj.json '. |= . + {"IntendedFor":$jsonSubjArr}' $fmapJson > ${derivDir}/tmpIntendAdd.json
				rm $fmapJson && mv ${derivDir}/tmpIntendAdd.json $fmapJson
			fi
		done
	fi
done
