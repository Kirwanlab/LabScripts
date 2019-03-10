#!/bin/bash


# Written by Nathan Muncy on 11/1/18


### --- Notes
#
# 1) This script will (a) read the summary files,
# 		(b) determine which subjects move too much
# 		(number of censored TRs >= 1.5*IQR  AND  > 0.03*TotalTRs)
# 		and (c) print a list of subjects to be removed.
#
#		R is needed for this script
#
# 2) Written for the supercomputer, but is light/fast
#
# 3) output written to $grpDir



module load r/3/5


###??? update these
parDir=~/compute/STT_reml
workDir=${parDir}/derivatives
grpDir=${parDir}/Analyses/grpAnalysis
refDir=${workDir}/sub-1295
maxP=0.1								# max percentage of censored TRs allowed (0.1 = 10%)

mkdir -p $grpDir



# make list of decons
cd $refDir

c=0; for i in *stats+tlrc.HEAD; do
	deconList[$c]=${i%_*}
	let c=$[$c+1]
done


### Determine subjects with many motion events
cd $workDir

for i in ${deconList[@]}; do

	# pull num of censored TRs
	print=${grpDir}/info_numCens_${i}.txt
	echo -e "Subj \t NumCen \t PropCen" > $print

	for j in s*; do

		cenN=`cat ${j}/out_summary_${i}.txt | grep -m 1 "TRs censored" | sed "s/[^0-9]//g"`
		cenP=`cat ${j}/out_summary_${i}.txt | grep "censor fraction" | awk '{print $4}'`
		echo -e "$j \t $cenN \t $cenP" >> $print
	done


	# determine max
	print1=${grpDir}/info_IQR_${i}.txt
	echo -e "IQR \t Max" > $print1

	echo `tail -n +2 $print | awk '{print $2}'` > tmp
	iqr=`Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'IQR(d)' < tmp | awk '{print $2}'`
	iqr15=`echo $iqr*1.5 | bc`

	quantiles=(`Rscript -e 'd<-scan("stdin", quiet=TRUE)' -e 'quantile(d)' < tmp`)
	max=`echo $iqr15+${quantiles[8]} | bc`
	echo -e "$iqr \t $max" >> $print1 && rm tmp


	# determine which subjs moved too much
	print2=${grpDir}/info_rmSubj_${i}.txt
	> $print2

	arrSubj=(`tail -n +2 $print | awk '{print $1}'`)
	arrNum=(`tail -n +2 $print | awk '{print $2}'`)
	arrPrp=(`tail -n +2 $print | awk '{print $3}'`)
	arrLen=${#arrSubj[@]}

	c=0; while [ $c -lt $arrLen ]; do
		if [ $(echo ${arrNum[$c]}'>'$max) == 1 ] || [ $(echo ${arrPrp[$c]}'>'$maxP | bc) == 1 ]; then

			echo ${arrSubj[$c]} >> $print2
		fi
		let c=$[$c+1]
	done
done
