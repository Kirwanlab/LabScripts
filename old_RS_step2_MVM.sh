#!/bin/bash


# written by Nathan Muncy on 12/12/17


# Notes:
#
# 1. This script starts with seed-based correlation matrices (FINAL_<seed>)
#       Z-transforms, constructs brain mask, runs MC simulations, and then an MVM for each seed.
#
# 2. Requires a reference dir to get functional voxel dimensions
#
# 3. Will create, and put output in $conDir
#
# 4. Updated for new ACF MC method (nate, 10/17/18)




workDir=/Volumes/K_Drive/MNI_RDoC/MRI_processed
tempDir=/Volumes/Yorick/Templates/vold2_mni
conDir=${workDir}/conAnalysis
refDir=${workDir}/P002  # Dir to reference voxel dimensions and files

mkdir -p $conDir


# Seed list, should match the number of seeds from previous script (JLF + coord seeds)
cd $refDir

c=0; for i in FINAL*HEAD; do

    tmp=${i##*_}
    roi=${tmp%+*}
    seedList[$c]=$roi

    let c=$[$c+1]
done



### Go
cd $workDir

print=${conDir}/acf_all.txt
> $print

for i in P*; do
cd $i

    ### Fischer's Z-transform
    for j in FINAL*HEAD; do

        scan=${j%.*}
        tmp=${j%+*}
        string=${tmp##*_}

        if [ ! -f ZTrans_${string}+tlrc.HEAD ]; then
            3dcalc -a $scan -expr 'log((1+a)/(1-a))/2' -prefix ZTrans_$string
        fi
    done


    ### pull noise calcs
    cat blur_errts.1D >> $print

cd $workDir
done




### make a binary brain mask - it would be better to use a gray matter mask
cd $conDir

if [ ! -f brain_mask_resampled+tlrc.HEAD ]; then

    cp ${tempDir}/priors_ACT/Template_BrainCerebellumBinaryMask.nii.gz ${conDir}/tmp_brain_mask.nii.gz
    3dcopy tmp_brain_mask.nii.gz tmp_brain_mask+tlrc
    3dfractionize -template ${refDir}/ZTrans_PCC+tlrc -prefix tmp_brain_mask_resampled -input tmp_brain_mask+tlrc
    3dcalc -a tmp_brain_mask_resampled+tlrc -prefix brain_mask_resampled -expr "step(a)" && rm tmp*
fi




### MC correction

# clean
sed '/ 0  0  0    0/d' acf_all.txt > acf_cleaned.txt

# calc mean
xA=`awk '{ total += $1 } END { print total/NR }' acf_cleaned.txt`
xB=`awk '{ total += $2 } END { print total/NR }' acf_cleaned.txt`
xC=`awk '{ total += $3 } END { print total/NR }' acf_cleaned.txt`

echo "$xA $xB $xC" > acf_average.txt

# simulate
mask=brain_mask_resampled+tlrc
3dClustSim -mask $mask -LOTS -iter 10000 -acf $xA $xB $xC > MCstats.txt




### Run MVMs, for each seed
# update 3dMVM for your own study

for a in ${seedList[@]}; do

    # continuous ISI
    scan=ZTrans_${a}+tlrc
    out=MVM_${a}

    if [ ! -f ${out}+tlrc.HEAD ]; then

        3dMVM -prefix $out -jobs 6 -mask brain_mask_resampled+tlrc \
        -bsVars 'ISI' \
        -qVars 'ISI' \
        -num_glt 1 \
        -gltLabel 1 ISI -gltCode 1 'ISI : ' \
        -dataTable \
        Subj    ISI     InputFile \
        P002	0       ${workDir}/P002/"$scan" \
        P003	6       ${workDir}/P003/"$scan" \
        P005	3       ${workDir}/P005/"$scan" \
        P007	12      ${workDir}/P007/"$scan" \
        P008	5       ${workDir}/P008/"$scan" \
        P009	8       ${workDir}/P009/"$scan" \
        P010	5       ${workDir}/P010/"$scan" \
        P011	7       ${workDir}/P011/"$scan" \
        P012	8       ${workDir}/P012/"$scan" \
        P013	1       ${workDir}/P013/"$scan" \
        P014	0       ${workDir}/P014/"$scan" \
        P015	3       ${workDir}/P015/"$scan" \
        P016	8       ${workDir}/P016/"$scan" \
        P019	1       ${workDir}/P019/"$scan" \
        P020	11      ${workDir}/P020/"$scan" \
        P021	4       ${workDir}/P021/"$scan" \
        P022	7       ${workDir}/P022/"$scan" \
        P025	6       ${workDir}/P025/"$scan" \
        P026	14      ${workDir}/P026/"$scan" \
        P027	4       ${workDir}/P027/"$scan" \
        P028	7       ${workDir}/P028/"$scan" \
        P029	14      ${workDir}/P029/"$scan" \
        P031	12      ${workDir}/P031/"$scan" \
        P033	19      ${workDir}/P033/"$scan" \
        P035	14      ${workDir}/P035/"$scan" \
        P037	1       ${workDir}/P037/"$scan"
    fi
done





