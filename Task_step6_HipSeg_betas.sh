#!/bin/bash





# Written by Nathan Muncy on 12/14/18


### --- Notes
#
# 1) This script will pull mean betas from L/R CA1, CA2/3/DG (called Multi), and Sub.
#		- maybe I'll update this in the future for other MTL regions
#		- also, maybe I'll update this to support more than 2 betas p/comparison
#				or...you could
#
# 2) Specifically, each mask for each hemisphere will be resampled,
#		binarized, and voxels defined by mutliple over-lapping masks
#		will be excluded.
#
# 3) Will also blur data (since we don't blur for ETAC)
#
# 4) A print out of the number of voxels in/excluded is supplied (info_*.txt)
#
# 5) Again, betas will not be extracted from participants who moved too much





# general vars											###??? Update these
workDir=/Volumes/Yorick/STT_new
roiDir=${workDir}/Analyses/roiAnalysis
grpDir=${workDir}/Analyses/grpAnalysis
priorDir=/Volumes/Yorick/Templates/vold2_mni/priors_HipSeg
refFile=${workDir}/s1295/SpT1_stats


# decon vars
compList=(SpT1 SpT1pT2 T1 T1pT2 T2 T2fT1)				# matches decon prefix
arrA=(33 39 53 59 97 103)								# setA beh sub-brik for etacList
arrB=(36 42 56 62 100 106)								# steB
compLen=${#compList[@]}
blurX=2													# blur multiplier




# function - search array for string
MatchString (){
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}




### make Sub, CA1, CA2/3/DG masks
mkdir -p $roiDir
cd $roiDir

# ref files
if [ ! -f ${refFile}.nii.gz ]; then
	3dcopy ${refFile}+tlrc ${refFile}.nii.gz
fi


gridSize=`fslhd ${refFile}.nii.gz | grep "pixdim1" | awk '{print $2}'`
int=`printf "%.0f" $gridSize`
blurInt=$(($blurX * $int))


if [ ! -f ${refFile}_blur${blurInt}+tlrc.HEAD ]; then
	3dmerge -prefix ${refFile}_blur${blurInt} -1blur_fwhm $blurInt -doall ${refFile}+tlrc
	3dcopy ${refFile}_blur${blurInt}+tlrc ${refFile}_blur${blurInt}.nii.gz
fi



for i in L R; do
	if [ ! -f Mask_${i}_CA1+tlrc.HEAD ]; then

		# resample
		for j in CA{1..3} DG Sub; do

			c3d ${priorDir}/${i}_${j}_prob.nii.gz -thresh 0.3 1 1 0 -o tmp_${i}_${j}.nii.gz
			3dfractionize -template ${refFile}_blur${blurInt}+tlrc -input tmp_${i}_${j}.nii.gz -prefix tmp_${i}_${j}_res
			3dcalc -a tmp_${i}_${j}_res+tlrc -prefix tmp_${i}_${j}_bin -expr "step(a-3000)"
			3dcopy tmp_${i}_${j}_bin+tlrc tmp_${i}_${j}_bin.nii.gz
		done


		# stitch 2/3/DG (Multi)
		c3d tmp_${i}_CA2_bin.nii.gz tmp_${i}_CA3_bin.nii.gz tmp_${i}_DG_bin.nii.gz -accum -add -endaccum -o tmp_${i}_Multi.nii.gz
		c3d tmp_${i}_Multi.nii.gz -thresh 0.1 inf 1 0 -o tmp_${i}_Multi_bin.nii.gz


		# exclude overlapping voxels
		c3d tmp_${i}_Multi_bin.nii.gz tmp_${i}_CA1_bin.nii.gz tmp_${i}_Sub_bin.nii.gz -accum -add -endaccum -o tmp_${i}_sum.nii.gz
		c3d tmp_${i}_sum.nii.gz -thresh 1.1 10 0 1 -o tmp_${i}_rm.nii.gz
		c3d tmp_${i}_sum.nii.gz -dup -lstat > info_${i}sum.txt

		for k in Multi Sub CA1; do
			c3d tmp_${i}_${k}_bin.nii.gz tmp_${i}_rm.nii.gz -multiply -o Mask_${i}_${k}.nii.gz
			c3d Mask_${i}_${k}.nii.gz -dup -lstat > info_${i}_${k}.txt
			3dcopy Mask_${i}_${k}.nii.gz Mask_${i}_${k}+tlrc
		done

		rm tmp* Mask*.nii.gz
	fi
done




### Pull Betas
for i in ${!compList[@]}; do

	pref=${compList[$i]}
	scan=${pref}_stats_blur${blurInt}+tlrc
	betas=${arrA[$i]},${arrB[$i]}

	arrRem=(`cat ${grpDir}/info_rmSubj_${pref}.txt`)
	print=${roiDir}/Betas_${pref}_sub_data.txt
	> $print

	for k in Mask*.HEAD; do

		hold=${k#*_}
		echo "Mask ${hold%+*}" >> $print

		for j in ${workDir}/s*; do

			subj=${j##*\/}
			MatchString $subj "${arrRem[@]}"

			if [ $? == 1 ]; then

				if [ ! -f ${j}/${scan}.HEAD ]; then
					3dmerge -prefix ${j}/${scan%+*} -1blur_fwhm ${blurInt} -doall ${j}/${scan%_*}+tlrc
				fi

				stats=`3dROIstats -mask ${k%.*} "${j}/${scan}[${betas}]"`
				echo "$subj $stats" >> $print
			fi
		done

		echo >> $print
	done
done


> Master_list.txt
for i in Betas*; do
	echo $i >> Master_list.txt
done
