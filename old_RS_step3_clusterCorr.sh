#!/bin/bash


# Written by Nathan Muncy on 12/18/17


# Notes:
#
# 1. Assumes all clusters (from MVM) are titled Clust_<seed> and exist in $conDir
#       Also, make sure that clusters are title according to seed names, in this format: Clust_<seedName>_mask i.e. Clust_rPCC_mask+tlrc
#
# 2. Will create new dirs in $conDir ($clustDir, $corrDir) and write output there
#
# 3. This will start with clusters, split them, pull participant correlations for each mask, and print them to $corrDir.
#       It will also pull cluster mask coordinates, and print them to $corrDir.



workDir=/Volumes/K_Drive/MNI_NREM/MRI_processed
conDir=${workDir}/conAnalysis
clustDir=${conDir}/mvm_clusters
betaDir=${conDir}/mvm_betas



# Organize
cd $conDir

mkdir $clustDir
mkdir $betaDir
mv Clust* $clustDir



# make array of seeds that had results
cd $clustDir

c=0; for i in Clust*HEAD; do
    if [[ $i == *mask+tlrc.HEAD ]]; then

        tmp=${i%_*}
        roi=${tmp##*_}
        arr[$c]=$roi

        let c=$[$c+1]
    fi
done



# Make subj list
cd $workDir

c=0; for i in P*; do
    subjList[$c]=$i
    let c=$[$c+1]
done




### Work
cd $clustDir

# make list of parent clusters
c=0; for i in Clust*mask+tlrc.HEAD; do
    clustList[$c]=${i%.*}
    let c=$[$c+1]
done


# Split clusters into masks, make list
for i in ${clustList[@]}; do

    numMasks=`3dinfo $i | grep "At sub-brick #0 '#0' datum type is short" | sed 's/[^0-9]*//g' | sed 's/^...//'`
    hold=${i%+*}
    3dcopy $i tmp_${hold}.nii.gz

    for (( c=1; c<=$numMasks; c++ )); do
        if [ ! -f ${hold}_${c}+tlrc.HEAD ]; then

            c3d tmp_${hold}.nii.gz -thresh $c $c 1 0 -o tmp_${hold}_${c}.nii.gz
            3dcopy tmp_${hold}_${c}.nii.gz ${hold}_${c}+tlrc
        fi
    done

rm tmp*
done



# Pull, print betas
for i in ${arr[@]}; do

    print=${betaDir}/Betas_${i}.txt
    > $print

    for j in Clust*_${i}_mask_*.HEAD; do
        for k in ${subjList[@]}; do

            subjDir=${workDir}/$k
            file=FINAL_Seed_${i}+tlrc
            mask=${j%.*}

            hold=${mask%+*}
            maskNum=mask${hold##*_}

            stats=`3dROIstats -mask $mask "${subjDir}/${file}[0]"`
            echo "$maskNum $k $stats" >> $print

        done

        echo >> $print
    done
done


# Pull mask coordinates
print1=${betaDir}/MaskCoordinates.txt
> $print1

cd $clustDir
for i in *table.1D; do

    tmp=${i%_*}
    seed=${tmp##*_}

    echo $seed >> $print1
    cat $i >> $print1
    echo >> $print1

done





























