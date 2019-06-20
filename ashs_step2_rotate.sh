#!/bin/bash


# Written by Nathan Muncy on 3/20/18




### Notes:
#
# 1. This script is written to rotate and resample ASHS masks into participant functional space
#
# 2. CA2, 3, and DG will be combined.
#
# 3. Specify which ashs output files you will use, for $Lashs and $Rashs



workDir=/Volumes/Yorick/<somewhere>

refFile=volreg+orig		     	# participant file with desired dim/rotation
Lashs=left_lfseg_corr_usegray  	# desired left ashs string
Rashs=right_lfseg_corr_usegray  # desired right ashs string

LB=0.3							# lower bound for thresholding (0-1)

numArr=(1 8 9 11 12)   			# desired ashs masks, in addition to CA23DG
namArr=(CA1 SUB ERC BA35 BA36)	# names for desired ashs masks



cd $workDir

for i in s*; do
cd $i

    if [ ! -f Lhipp_CA23DG.nii.gz ]; then

		# get data
        cp final/${i}_${Lashs}.nii.gz ./tmp_Lhipp.nii.gz
        cp final/${i}_${Rashs}.nii.gz ./tmp_Rhipp.nii.gz


        # split labels to avoid blurring, rotate, resample, thresh
        for j in L R; do
            for k in {1..13}; do

                c3d tmp_${j}hipp.nii.gz -thresh $k $k 1 0 -o tmp_${j}hipp_${k}.nii.gz
                3dwarp -oblique_parent $refFile -prefix tmp_${j}hipp_${k}_rotated.nii.gz -gridset $refFile -quintic tmp_${j}hipp_${k}.nii.gz
                c3d tmp_${j}hipp_${k}_rotated.nii.gz -thresh $LB 1 1 0 -o ${j}hipp_${k}.nii.gz
            done
        done


        # combine CA2,3,DG and binarize
        for a in L R; do

            c3d ${a}hipp_2.nii.gz ${a}hipp_3.nii.gz ${a}hipp_4.nii.gz -accum -add -endaccum -o tmp_${a}hipp_CA23DG.nii.gz
            c3d tmp_${a}hipp_CA23DG.nii.gz -thresh 0.1 10 1 0 -o ${a}hipp_CA23DG.nii.gz
        done


        # rename output files
        for a in Lhipp* Rhipp*; do
            num="${a//[^0-9]/}"

            for b in ${!numArr[@]}; do
                if [[ $num == ${numArr[$b]} ]]; then

                    name=${namArr[$b]}
                    mv $a "${a/$num/$name}"
                fi
            done
        done


		# clean up unused
        rm tmp*
        for x in _{2,3,4,5,6,7,10,13}; do
            rm *${x}*
        done
    fi

cd $workDir
done










