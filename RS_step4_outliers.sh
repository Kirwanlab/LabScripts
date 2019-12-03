#!/bin/bash


# written by Nathan Muncy on 12/3/19


### --- Notes
#
# 1) This script will (a) read the summary files,
# 	(b) determine which subjects move too much
# 	(number of censored TRs >= 1.5*IQR  AND/OR  > 0.03*TotalTRs)
# 	and (c) print a list of subjects to be removed.
#
#	R is needed for this script
#
# 2) Written for the supercomputer, but is light/fast
#
# 3) Output written to $grpDir
#
# 4) Assumes only 1 out_summary.txt file exists for e/participant
#
# Researcher Input A) Array of phases/session
# Researcher Input B) max percentage of censored volumes allowed (0.1 = 10%)




### --- Set Up --- ###

parDir=~/compute/MNI_NREM/MRI_processed   	# working directory
workDir=${parDir}/derivatives				# derivative directory
grpDir=${parDir}/Analyses/grpAnalysis		# output directory for group-level analysis


phaseArr=(Sleep Wake)						###??? Researcher Input A
maxP=0.1									###??? Researcher Input B




### --- Make exclusion List --- ###
#
# Determine which participant moved excessively.
#
# IQR reported in info_$phase_IQR.txt
# Table of movement in info_$phase_numCens.txt
# List of subjects to be excluded in info_$phase_rmSubj.txt


mkdir -p $grpDir
module load r/3.6


cd $workDir

for i in ${phaseArr[@]}; do

	# pull num of censored TRs
	print=${grpDir}/info_${i}_numCens.txt
	echo -e "Subj \t NumCen \t PropCen" > $print

	for j in s*; do

		cenN=`cat ${j}/ses-${i}/out_summary.txt | grep -m 1 "TRs censored" | sed "s/[^0-9]//g"`
		cenP=`cat ${j}/ses-${i}/out_summary.txt | grep "censor fraction" | awk '{print $4}'`
		echo -e "$j \t $cenN \t $cenP" >> $print
	done


	# determine max according to IQR method
	print1=${grpDir}/info_${i}_IQR.txt
	echo -e "IQR \t Max" > $print1

	echo `tail -n +2 $print | awk '{print $2}'` > tmp
	iqr=`Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'IQR(d)' < tmp | awk '{print $2}'`
	iqr15=`echo $iqr*1.5 | bc`

	quantiles=(`Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'quantile(d)' < tmp`)
	max=`echo $iqr15+${quantiles[8]} | bc`
	echo -e "$iqr \t $max" >> $print1 && rm tmp


	# determine which subjs moved too much
	print2=${grpDir}/info_${i}_rmSubj.txt
	> $print2

	arrSubj=(`tail -n +2 $print | awk '{print $1}'`)
	arrNum=(`tail -n +2 $print | awk '{print $2}'`)
	arrPrp=(`tail -n +2 $print | awk '{print $3}'`)
	arrLen=${#arrSubj[@]}

	c=0; while [ $c -lt $arrLen ]; do
		# if [ $(echo ${arrNum[$c]}'>'$max | bc) == 1 ] && [ $(echo ${arrPrp[$c]}'>'$maxP | bc) == 1 ]; then
		if [ $(echo ${arrPrp[$c]}'>'$maxP | bc) == 1 ]; then

			echo ${arrSubj[$c]} >> $print2
		fi
		let c=$[$c+1]
	done
done

