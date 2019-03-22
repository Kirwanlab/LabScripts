#!/bin/bash


# Written by Nathan Muncy on 3/20/18


# Notes:
#
# 1. This script will read the text file of coordinates, and create a single-voxel
#		mask at the peak value for each mask. It will then extract the correlation
#		of each subject at that coordinate.
#
# 2. Manual intervention is needed at the "Split Coordinates" step, where you strip
#		off the top and bottom of the text file, around the useless parts. This could
#		be supplemented in R.






### Set up
workDir=/Volumes/K_Drive/MNI_NREM/MRI_processed
conDir=${workDir}/conAnalysis
betaDir=${conDir}/mvm_betas
clustDir=${conDir}/mvm_clusters
peakDir=${conDir}/mvm_peak


if [ ! -d $peakDir ]; then
    mkdir $peakDir
fi




### Function
MakeVox () {

    name=$1

    eval arrX=("$2")
    eval arrY=("$3")
    eval arrZ=("$4")
    arrLen=${#arrX[@]}

    cc=0; dd=1; while [ $cc -lt $arrLen ]; do
        echo "${arrX[$cc]} ${arrY[$cc]} ${arrZ[$cc]}" > tmp.txt
        3dUndump -prefix ${name}_mask${dd} -master ${conDir}/MVM_CvI_${name}+tlrc -xyz tmp.txt
        rm tmp.txt
        cc=$dd; dd=$[$dd+1]
    done
}




### Master array, of all seeds
cd $clustDir

c=0; for i in Clust*HEAD; do
    if [[ $i == *mask+tlrc.HEAD ]]; then

        tmp=${i%_*}
        roi=${tmp##*_}
        seedList[$c]=$roi

        let c=$[$c+1]
    fi
done

seedLen=${#seedList[@]}

## Check syntax
# echo ${seedList[@]}



### Split Coordinates, for all seeds in MaskCoordinates.txt
# --- CHANGE THIS FOR YOUR STUDY --- #
cd $betaDir

if [ ! -f coord_lLing.txt ]; then

	sed '1,6d' MaskCoordinates.txt | sed '4,$d' > coord_lLing.txt
	sed '1,16d' MaskCoordinates.txt | sed '2,$d' > coord_lMFG.txt
	sed '1,24d' MaskCoordinates.txt | sed '5,$d' > coord_lPar.txt
	sed '1,35d' MaskCoordinates.txt | sed '6,$d' > coord_rLing.txt
	sed '1,47d' MaskCoordinates.txt > coord_rPCC.txt
fi



### Make peak masks (1 voxel) for each cluster
cd $peakDir

c=0; while [ $c -lt $seedLen ]; do

    string=${seedList[$c]}

    # Coord array for each seed
    if [ ! -f ${string}_mask1+tlrc.HEAD ]; then

		holdX=($(cat ${betaDir}/coord_${string}.txt | awk '{print $5}'))
		holdY=($(cat ${betaDir}/coord_${string}.txt | awk '{print $6}'))
		holdZ=($(cat ${betaDir}/coord_${string}.txt | awk '{print $7}'))

		# make single voxel mask for e/set of coord
		MakeVox $string "$(echo ${holdX[@]})" "$(echo ${holdY[@]})" "$(echo ${holdZ[@]})"
	fi


    # Fill master array (seed[0]=seed_mask1+tlrc), with output of MakeVox
    d=0; for i in ${string}_mask*.HEAD; do

        eval ${string}[$d]=${i%.*}

    let d=$[$d+1]
    done

let c=$[$c+1]
done

## Check syntax
# for i in ${!seedList[@]}; do
# 	eval echo \${${seedList[$i]}[@]}
# done



### Participant peak R
cd $workDir

# Pull nested-array data; extract voxel beta from appropriate file for each mask
for i in ${!seedList[@]}; do

	seed=${seedList[$i]}
	scan=FINAL_Seed_${seed}+tlrc

	print=${peakDir}/${seed}_peak.txt
	> $print

	eval maskArr=(\${${seedList[$i]}[@]})

	## Check syntax
# 	echo ${maskArr[@]}

	for j in ${maskArr[@]}; do

		mask=${peakDir}/$j

		for k in P*; do     # loop through participants in $workDir

			file=${workDir}/${k}/$scan
			stats=`3dROIstats -mask $mask "${file}[0]"`
			echo "$k $j $stats" >> $print
		done

		echo >> $print
	done
done




