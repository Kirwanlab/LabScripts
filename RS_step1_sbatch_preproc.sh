#!/bin/bash

#SBATCH --time=90:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16gb   # memory per CPU core
#SBATCH -J "RSF"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE





# This is derived from afni_proc.py: Example 10
# Written by Nathan Muncy on 11/7/17
# Bandpass regressor file creation and inclusion in regression added by Ty Bodily on 1/23/18

# Notes:
# 1. This assumes that functional files exist in NIfTI format (Resting.nii.gz), uses this header data
#       to determine TR count and voxel size
#
# 2. Uses AFNI, FSL, c3d, and ANTs commands
#
# 3. Is built to render ROI (JLF) or coordinate seeds
#
# 4. Assumes that the relevant JLF/ACT priors exist in the template dir
#
# 5. This script starts with aw 3D struct and 4D RS files, and finishes with seed-based correlation matrices (FINAL_*, Pearson's R).
#
# 6. For seeds, the script will name each input (e.g. roiName corresponds to roiList). These two arrays must have an equal length.



### Set variables
# General variables
tempDir=~/bin/Templates/vold2_mni    # Location of template/priors
actDir=${tempDir}/priors_ACT
jlfDir=${tempDir}/priors_JLF
scriptDir=~/compute/Scripts    # Location of antify script


# Subject variables
workDir=${1}/$2


# Priors for WM/CSF masks
priorList=(3 0015 0014 0004 0043)
labelList=(WM 4ven 3ven LLven RLven)
arrLen=${#labelList[@]}



### Study-specific variables
# JLF seeds
roiList=(2025)
roiName=(jlf_PCC)
roiLen=${#roiList[@]}

# Coordinate seeds
coordList=("5 -55 25" "-40 22 38" "-30 -52 54" "-26 -72 -14" "8 -78 -6")
coordName=(rPCC lMFG lPar lLing rLing)
coordLen=${#coordList[@]}
seedSize=5      # Size (mm) of radius

# Flags
jlfFlag=0
corFlag=0
blurFlag=0




### Make it happen
cd $workDir

# Get outlier info
if [ ! -f outcount_Resting.1D ]; then

    3dcopy Resting.nii.gz Resting+orig
    3dToutcount -automask -fraction -polort 5 -legendre Resting+orig > outcount_Resting.1D
    1deval -a outcount_Resting.1D -expr "1-step(a-0.1)" > outcount_censor.1D
fi



# Get var info
tr_count=`fslhd Resting.nii.gz | grep "dim4" | awk 'FNR == 1 {print $2}'`
minindex=`3dTstat -argmin -prefix - outcount_Resting.1D\'`
ovals=(`1d_tool.py -set_run_lengths $tr_count -index_to_run_tr $minindex`)

minoutrun=${ovals[0]}
minouttr=${ovals[1]}



# Despike, get min outlier volume, volreg
if [ ! -f volreg_Resting+orig.HEAD ]; then

    3dDespike -NEW -nomask -prefix despike_Resting Resting+orig
    3dbucket -prefix min_outlier_volume despike_Resting+orig"[$minouttr]"     # new base for volreg
    3dvolreg -zpad 1 -base min_outlier_volume+orig -1Dfile dfile_Resting.1D -prefix volreg_Resting -cubic despike_Resting+orig
fi



# Construct motion files
if [ ! -f censor_combined_2.1D ]; then

    1d_tool.py -infile dfile_Resting.1D -set_nruns 1 -demean -write motion_demean.1D
    1d_tool.py -infile dfile_Resting.1D -set_nruns 1 -derivative -demean -write motion_deriv.1D
    1d_tool.py -infile dfile_Resting.1D -set_nruns 1 -show_censor_count -censor_prev_TR -censor_motion 0.2 motion  # more conservative threshold should be used for RSfunc
    1deval -a motion_censor.1D -b outcount_censor.1D -expr "a*b" > censor_combined_2.1D
fi



# ANTs reg
dim=3
out=ants_
fix=${tempDir}/vold2_mni_head.nii.gz
mov=struct_rotated.nii.gz

if [ ! -f ants_0GenericAffine.mat ]; then

    3dWarp -oblique_parent volreg_Resting+orig -prefix $mov struct_raw.nii.gz

    antsRegistrationSyN.sh \
    -d $dim \
    -f $fix \
    -m $mov \
    -o $out
fi



# Antify
if [ ! -f volreg_Resting_ANTS_resampled+tlrc.HEAD ]; then
    ${scriptDir}/antifyFunctional_nate.sh -a $out -t $fix -i volreg_Resting+orig
fi



# Get masks, resample
if [ ! -f mask_WM+tlrc.HEAD ]; then
    c=0; while [ $c -lt $arrLen ]; do

        if [ ${priorList[$c]} == 3 ]; then
            prior=${actDir}/Prior${priorList[$c]}.nii.gz
        else
            prior=${jlfDir}/label_${priorList[$c]}.nii.gz
        fi

        name=mask_${labelList[$c]}.nii.gz
        3dWarp -oblique_parent volreg_Resting_ANTS_resampled+tlrc -prefix tmp_$name -gridset volreg_Resting_ANTS_resampled+tlrc -quintic $prior

    let c=$[$c+1]
    done


    c3d tmp*ven.nii.gz -accum -add -endaccum -o tmp_mask_CSF.nii.gz

    for a in tmp_mask_{WM,CSF}; do

        3dcopy ${a}.nii.gz ${a}+tlrc
        3dcalc -a ${a}+tlrc -prefix ${a#*_}+tlrc -expr "ispositive(a-0.5)"

    done
    rm tmp*
fi



# Tissue-class TimeSeries
if [ ! -f tissClass_WM.1D ]; then
    for j in CSF WM; do
        3dmaskave -quiet -mask mask_${j}+tlrc volreg_Resting_ANTS_resampled+tlrc | 1d_tool.py -infile - -demean -write tissClass_${j}.1D
    done
fi



# Note non-censored TRs
ktrs=`1d_tool.py -infile censor_combined_2.1D -show_trs_uncensored encoded`



# get exclusion mask
if [ ! -f full_mask+tlrc.HEAD ]; then

    cp ${actDir}/Template_BrainCerebellumBinaryMask.nii.gz ./brain_mask.nii.gz
    3dfractionize -template volreg_Resting_ANTS_resampled+tlrc -input brain_mask.nii.gz -prefix full_mask
fi



#create bandpass regressor 1d file. Period from .008 to .08 Hz spans 12.5-125 s.
1dBport -nodata $tr_count -band 0.008 0.08 -nozero -invert > bandpass_rall.1D




# Deconvolve - remove noise from data
if [ ! -f errts+tlrc.HEAD ]; then

    3dDeconvolve -input volreg_Resting_ANTS_resampled+tlrc \
    -mask full_mask+tlrc \
    -censor censor_combined_2.1D \
    -ortvec tissClass_WM.1D tiss.WM \
    -ortvec tissClass_CSF.1D tiss.CSF \
    -ortvec bandpass_rall.1D bandpass \
    -polort A \
    -num_stimts 12 \
    -stim_file 1  motion_demean.1D'[0]' -stim_base 1  -stim_label 1 roll_01  \
    -stim_file 2  motion_demean.1D'[1]' -stim_base 2  -stim_label 2 pitch_01 \
    -stim_file 3  motion_demean.1D'[2]' -stim_base 3  -stim_label 3 yaw_01   \
    -stim_file 4  motion_demean.1D'[3]' -stim_base 4  -stim_label 4 dS_01    \
    -stim_file 5  motion_demean.1D'[4]' -stim_base 5  -stim_label 5 dL_01    \
    -stim_file 6  motion_demean.1D'[5]' -stim_base 6  -stim_label 6 dP_01    \
    -stim_file 7  motion_deriv.1D'[0]'  -stim_base 7  -stim_label 7 roll_02  \
    -stim_file 8  motion_deriv.1D'[1]'  -stim_base 8  -stim_label 8 pitch_02 \
    -stim_file 9  motion_deriv.1D'[2]'  -stim_base 9  -stim_label 9 yaw_02   \
    -stim_file 10 motion_deriv.1D'[3]'  -stim_base 10 -stim_label 10  dS_02  \
    -stim_file 11 motion_deriv.1D'[4]'  -stim_base 11 -stim_label 11  dL_02  \
    -stim_file 12 motion_deriv.1D'[5]'  -stim_base 12 -stim_label 12  dP_02  \
    -fout -tout -x1D X.xmat.1D -x1D_uncensored X.nocensor.xmat.1D \
    -fitts fitts \
    -errts errts \
    -jobs 6
fi



# Blur
# done after deconvolution, so WM signal is not mixed with GM

hold=`fslhd Resting.nii.gz | grep "pixdim1" | awk '{print $2}'`
int=`printf "%.0f" $hold`
vsz="$(($int * 2))"
blurName=errts_out+tlrc

if [ $blurFlag == 1 ]; then
    if [ ! -f ${blurName}.HEAD ]; then

        3dmerge -1blur_fwhm $vsz -doall -prefix ${blurName%+*} errts+tlrc
    fi
else
    cp errts+tlrc.HEAD ${blurName}.HEAD
    cp errts+tlrc.BRIK ${blurName}.BRIK
fi


# Project out regression matrix
if [ ! -f errts_tproject+tlrc.HEAD ]; then

    3dTproject -polort 0 -input volreg_Resting_ANTS_resampled+tlrc -censor censor_combined_2.1D -cenmode ZERO \
    -ort X.nocensor.xmat.1D -prefix errts_tproject
fi



#1d_tool.py -show_cormat_warnings -infile X.xmat.1D |& tee out.cormat_warn.txt



# Combine runs, get signal to noise ratio
if [ ! -f TSNR+tlrc.HEAD ]; then

    3dTcat -prefix all_runs volreg_Resting_ANTS_resampled+tlrc
    3dTstat -mean -prefix rm.signal.all all_runs+tlrc"[$ktrs]"
    3dTstat -stdev -prefix rm.noise.all errts_tproject+tlrc"[$ktrs]"
    3dcalc -a rm.signal.all+tlrc -b rm.noise.all+tlrc -c full_mask+tlrc -expr 'c*a/b' -prefix TSNR
fi



# compute global correlation average
if [ ! -f out_gcor.1D ]; then

    3dTnorm -norm2 -prefix rm.errts.unit errts_tproject+tlrc
    3dmaskave -quiet -mask full_mask+tlrc rm.errts.unit+tlrc > gmean_errts_unit.1D
    3dTstat -sos -prefix - gmean_errts_unit.1D\' > out_gcor.1D
fi



# Compute correlation volume
if [ ! -f rm.DP+tlrc.HEAD ]; then

    3dcalc -a rm.errts.unit+tlrc -b gmean_errts_unit.1D -expr 'a*b' -prefix rm.DP
    3dTstat -sum -prefix corr_brain rm.DP+tlrc
fi



# Compute sum of non-baseline regressors
if [ ! -f sum_ideal.1D ]; then

    reg_cols=`1d_tool.py -infile X.nocensor.xmat.1D -show_indices_interest`
    3dTstat -sum -prefix sum_ideal.1D X.nocensor.xmat.1D"[$reg_cols]"
    1dcat X.nocensor.xmat.1d"[$reg_cols]" > X.stim.xmat.1D
fi


### 3dFWHMx
> blur_est.1D
> blur_epits.1D
> blur_errts.1D

# model uncensored TRs
trs=`1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded -show_trs_run 1`
3dFWHMx -detrend -mask full_mask+tlrc all_runs+tlrc"[${trs}]" >> blur_epits.1D
3dFWHMx -detrend -mask full_mask+tlrc errts_tproject+tlrc"[${trs}]" >> blur_errts.1D

# compute avg blur, append
blursA=`cat blur_epits.1D`
blursB=`cat blur_errts.1D`

echo "$blursA  # epits blur estimates" >> blur_est.1D
echo "$blursB  # errts blur estimages" >> blur_est.1D


# pull JLF label connectivity
if [ $jlfFlag == 1 ]; then
    c=0; while [ $c -lt $roiLen ]; do

        seedName=Seed_${roiName[$c]}

        if [ ! -f FINAL_${seedName}+tlrc.HEAD ]; then

            jlfPrior=${jlfDir}/label_${roiList[$c]}.nii.gz

            3dcopy $jlfPrior tmp_${seedName}+tlrc
            3dfractionize -template $blurName -input tmp_${seedName}+tlrc -prefix tmp_rs_${seedName}
            3dcalc -a tmp_rs_${seedName}+tlrc -prefix $seedName -expr "step(a)" && rm tmp*

            3dROIstats -quiet -mask ${seedName}+tlrc $blurName > ${seedName}_TimeSeries.1D
            3dTcorr1D -mask full_mask+tlrc -prefix FINAL_${seedName} $blurName ${seedName}_TimeSeries.1D
        fi
    let c=$[$c+1]
    done
fi



# generate ROI from coordinates, pull connectivity
if [ $corFlag == 1 ]; then
    c=0; while [ $c -lt $coordLen ]; do

        seedName=Seed_${coordName[$c]}

        if [ ! -f FINAL_${seedName}+tlrc.HEAD ]; then

            echo ${coordList[$c]} > ${seedName}.txt
            3dUndump -prefix $seedName -master $blurName -srad $seedSize -xyz ${seedName}.txt
            3dROIstats -quiet -mask ${seedName}+tlrc $blurName > ${seedName}_TimeSeries.1D
            3dTcorr1D -mask full_mask+tlrc -prefix FINAL_${seedName} $blurName ${seedName}_TimeSeries.1D
        fi
    let c=$[$c+1]
    done
fi












