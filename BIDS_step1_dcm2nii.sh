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




###??? change these variables/arrays
rawDir=/Volumes/Yorick/MriRawData							# location of raw data
workDir=/Volumes/Yorick/STT_reml							# desired working directory
tempDir=/Volumes/Yorick/Templates

subjList=(1295 1859 1875 1949 2439 2575 3013 3018 3023 3031 3042 3043 3046 3075 3223 3425 3710 3711 3712 3735 3746 3761 3762 3763 3767 3768 3775 3778 3783 3786 3816 3849 3850 3851 3852)									# list of subjects
session=STT													# scanning session - for raw data organization (ses-STT)
task=task-stt												# name of task, for epi data naming

epiDirs=(STUDY TEST1{1,2} TEST2{1..4})						# epi dicom directory name/prefix
t1Dir=t1													# t1 ditto
t2Dir=HHR													# t2 ditto
blipDir=()													# blip ditto - names one per scanning phase




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
	echo -e "\te.g. \"Nate Muncy\", \"Brock Kirwan\"\n"
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
	anatDir=${workDir}/rawdata/sub-${i}/anat
	funcDir=${workDir}/rawdata/sub-${i}/func
	derivDir=${workDir}/derivatives/sub-${i}

	for j in {anat,func,deriv}Dir; do
		if [ ! -d $(eval echo \${$j}) ]; then
			mkdir -p $(eval echo \${$j})
		fi
	done


	### construct data
	dataDir=${rawDir}/sub-${i}/ses-${session}/dicom

	# t1 data
	if [ ! -f ${anatDir}/sub-${i}_T1w.nii.gz ]; then

		dcm2niix -b y -ba y -z y -o $anatDir -f tmp_sub-${i}_T1w ${dataDir}/${t1Dir}*/

		# Deface by Brock Kirwan
		3dAllineate -base ${anatDir}/tmp_sub-${i}_T1w.nii.gz -input ${tempDir}/mean_reg2mean.nii.gz -prefix ${anatDir}/tmp_mean_reg2mean_aligned.nii -1Dmatrix_save ${anatDir}/tmp_allineate_matrix
		3dAllineate -base ${anatDir}/tmp_sub-${i}_T1w.nii.gz -input ${tempDir}/facemask.nii.gz -prefix ${anatDir}/tmp_facemask_aligned.nii -1Dmatrix_apply ${anatDir}/tmp_allineate_matrix.aff12.1D
		3dcalc -a ${anatDir}/tmp_facemask_aligned.nii -b ${anatDir}/tmp_sub-${i}_T1w.nii.gz -prefix ${anatDir}/sub-${i}_T1w.nii.gz -expr "step(a)*b"
		mv ${anatDir}/tmp_sub-${i}_T1w.json ${anatDir}/sub-${i}_T1w.json
		rm ${anatDir}/tmp*
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

		# Json append by Brock Kirwan
		funcJson=${funcDir}/sub-${i}_${task}_run-${pos}_bold.json
		taskExist=$(cat $funcJson | jq '.TaskName')
		if [ "$taskExist" == "null" ]; then
			jq '. |= . + {"TaskName":"STT"}' $funcJson > ${derivDir}/tasknameadd.json         ###??? replace "STT" with your task name
			rm $funcJson && mv ${derivDir}/tasknameadd.json $funcJson
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
