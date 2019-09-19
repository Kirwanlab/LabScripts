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
# at some point for a variable duration.
#
# 4) TimingNames is only used for 1D files. List these behaviors in the same order that bash would list the timing files. E.g.:
#
#		TimingNames=(Hit Miss) would correspond to
#		behVect.01.1D (Hit) and behVect.02.1D (Miss)
#
# 5) Naming for txt files should be embedded in the file name. E.g. sub-1234_TF_Hit_All.txt
#		Determining behavior ($beh) for txt files is currently written for Temporal timing files named like the previous example
#
# 6) Assumes that the only files in derivatives/sub-1234/timing_files are the 1D/txt timing files



###??? update here
parDir=~/compute/Temporal/Experiment3			# research directory
derivDir=${parDir}/derivatives
rawDir=${parDir}/rawdata


taskName=Temporal								# name of task (ref Task_step0)
duration=										# duration of trial (ref Task_step2). Leave empty of duration varies (e.g. 1:1.3 5:0.8)
suffix=txt		   								# suffix of timing file (1D or txt)
TimingNames=(CR Hit LureFA LureCR Other) 		# behaviors corresponding to 1D timing files (ref Task_step2 NotLazy section)




### --- Work --- ###
#
# Makes sure that the correct number of timing files and runs exist.
# Writes a sorted timing file for each run in BIDs format


cd $rawDir

for i in sub-*; do

	# Set subj paths
	tfPath=${derivDir}/${i}/timing_files
	rawPath=${rawDir}/${i}/func


	# Determine number of runs, timing files, timing file rows
	epiCount=`ls ${rawPath}/*task*.nii.gz | wc -l`
	tfArr=(`ls ${tfPath}/*.$suffix | sed 's/.*\///'`)
	tfRow=`cat ${tfPath}/${tfArr[0]} | wc -l`

	# Check
	if [ $epiCount != $tfRow ]; then
		echo "Unequal number of EPI runs and rows of timing files for $i. Exit 1" >&2
		exit 1
	fi


	## extract rows from each TF, add beh for each run/$tfRow
	for((j=1; j<=$tfRow; j++)); do

		# start tmp file for each run/row
		tmpFile=${rawPath}/tmp_event_run-${j}.txt
		> $tmpFile

		## extract row from each TF, write it to $tmpFile & add beh, duration
		c=0; while [ $c -lt ${#tfArr[@]} ]; do

			# determine behavior from $TimingNames or txt file name
			if [ $suffix == 1D ]; then
				beh=${TimingNames[$c]}
			elif [ $suffix == txt ]; then
				tmp=${tfArr[$c]#*TF_}; beh=${tmp%_*}								### This is for Temporal timing files
			fi

			# pull values from row $j of timing file
			holdArr=(`sed "${j}q;d" ${tfPath}/${tfArr[$c]}`)
			for k in ${holdArr[@]}; do

				# determine start, duration
				if [ -z $duration ]; then
					start=${k%\:*}
					dur=${k#*\:}
				else
					start=$k
				fi

				# write columns
				echo -e "${start}\t${dur}\t$beh" >> $tmpFile
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

