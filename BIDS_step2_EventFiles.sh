#!/bin/bash


# Written by Nathan Muncy on 9/7/19

# --- Notes
#
# 1) This script is written to run following a deconvolution script (e.g. Task_step2).
#
# 2) It will find timing files in derivatives/sub-123/timing_files, and use these
# to write out event files in rawdata/sub-123/func.
#
# 3) It is currently only written for fixed durations, an update will follow 
#at some point for a variable duration.




###??? update here
parDir=/Volumes/Yorick/Exercise
derivDir=${parDir}/derivatives
rawDir=${parDir}/rawdata


tf1d=1		   								# toggle for whether timing files are .1D (1) or .txt (0)
TFstring=behvect_correct 					# prefix string for timing files
tfNames=(CR Hit LureFA LureCR Other) 		# behaviors corresponding to timing files (ref Task_step2 NotLazy section)
taskName=mst								# name of task (ref Task_step0)
duration=3.5								# duration of trial (ref Task_step2)




### --- Work --- ###
#
# Makes sure that the correct number of timing files and runs exist.
# Writes a sorted timing file for each run in BIDs format


# determine TF suffix
if [ $tf1d == 1 ]; then
	suffix=1D
else
	suffix=txt
fi


# work
cd $rawDir

for i in s*; do

	# Set subj paths
	tfPath=${derivDir}/${i}/timing_files
	rawPath=${rawDir}/${i}/func


	# Determine number of runs, timing files, timing file rows
	epiCount=`ls ${rawPath}/*task*.nii.gz | wc -l`

	tfArr=(`ls ${tfPath}/${TFstring}*.$suffix | sed 's/.*\///'`)
	tfRow=`cat ${tfPath}/${tfArr[0]} | wc -l`


	# Check
	if [ $epiCount != $tfRow ]; then
		echo "Unequal number of EPI runs and rows of timing files for $i. Exit 1" >&2
		exit 1
	fi

	if [ ${#tfNames[@]} != ${#tfArr[@]} ]; then
		echo "Unequal number of Behaviors (\$tfNames) and Timing Files for $i. Exit 2" >&2
		exit 2
	fi


	## extract rows from each TF, add beh for each run/$tfRow
	for((j=1; j<=$tfRow; j++)); do

		# start tmp file for each run/row
		tmpFile=${rawPath}/tmp_event_run-${j}.txt
		> $tmpFile

		# extract row from each TF, write it to $tmpFile & add beh, duration
		c=0; while [ $c -lt ${#tfArr[@]} ]; do

			beh=${tfNames[$c]}
			holdArr=(`sed "${j}q;d" ${tfPath}/${tfArr[$c]}`)

			for k in ${holdArr[@]}; do
				echo -e "${k}\t${duration}\t$beh" >> $tmpFile
			done

			let c=$[$c+1]
		done


		# make header for event file
		eventFile=${rawPath}/${i}_task-${taskName}_run-${j}_events.tsv
		echo -e "onset\tduration\ttrial_type" > $eventFile

		# sort tmp file, write event file
		sort -k1 -n $tmpFile >> $eventFile
		rm $tmpFile
	done
done

