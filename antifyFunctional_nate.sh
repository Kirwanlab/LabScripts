#!/bin/bash



#################################
#####  Notes and Usage       ####
#################################


# antifyFunctional_nate.sh
#
# Script to apply ANTS transformation to afni functional data.
#
# Original perl script by Craig Stark
# Modified by Brock Kirwan
# Modified by Nathan Muncy
#
# Update 10/12/2015: Added affine transformation to the WarpImageMultiTransform call
#                    Added "deoblique" step before the afni to nii conversion to deal with alignment issues
#        11/19/2016: Ported to bash
#        01/10/2017: Patched for compatibility with antsRegistrationSyN.sh on MLF superComp
#                    Patched for compatibilty with the newer version of ANI:
#                        Precompiled binary linux_openmp_64: Nov 23 2016 (Version AFNI_16.3.13)
#                    Updated help
#                    Updated arguments
#                    Wrote detailed printout for MLF superComp
#                    Edited step 6 to get positive voxel dimensions





### Usage
function Usage {
cat << USAGE

##############################

 USAGE:  `basename $0` -a deformation_matrix -t model -i functional_dataset -o out_path -f out_file

 antifyFunctional warps functional data into template space, using calculations previously performed
 by registering T1-structural files to a template. It requires both the Warp.nii and Affine.txt
 from ANTs registration. NOTE: this version now calls for prefix_1Warp.nii.gz and prefix_0GenericAffine.mat,
 produced by antsRegistrationSyN.sh.

 Required Arguments:

    -a:     deformation_matrix: the prefix string from antsRegistrationSyN.sh (e.g. ants_ ).

    -t:     model: the absolute path to your model template.

    -i:     functional_dataset: the prefix string of the input functional dataset.


 Optional Arguments:

    -o:     out_path: the absolute path to the directory where a txt file will print, default = current working dir

    -f:     out_file: the name of the desired output txt file (e.g. s1234_output.txt), default = antify_output.txt


 EXAMPLE: `basename $0` -a ants_ -t /Volumes/Yorick/Templates/template.nii.gz -i func_data+orig -o /Volumes/Yorick/Study -f output.txt

##############################

USAGE
    exit 1
}




### Set Path
function SetPath {
cat << SETPATH

--------------------------------------------------------------------------------------
Error locating ANTS
--------------------------------------------------------------------------------------
It seems that the ANTSPATH environment variable is not set. Please add the ANTSPATH
variable. This can be achieved by editing the .bash_profile in the home directory.
Add:

ANTSPATH=/home/yourname/bin/ants/

Or the correct location of the ANTS binaries.

SETPATH
    exit 1
}




### Script Check
function BinCheck {
cat << BINCHECK

--------------------------------------------------------------------------------------
Error locating antsRegistrationSyN.sh
--------------------------------------------------------------------------------------

This script is written to be used with antsRegistrationSyN.sh. Please make sure this
script is in your ANTSPATH, and that you call for it in your script.

BINCHECK
    exit 1
}




#################################
#####  Set up                ####
#################################

# required args
deform_matrix=0
model=0
func_data=0

# opt args
out_path=$(pwd)
out_file=antify_output.txt

while getopts "a:f:i:o:t:" OPT; do
    case $OPT in

        a)
        deform_matrix=$OPTARG
        ;;

        f)
        out_file=$OPTARG
        ;;

        i)
        func_data=$OPTARG
        ;;

        o)
        out_path=$OPTARG
        ;;

        t)
        model=$OPTARG
        ;;

    esac
done


tmp=`echo $func_data | sed 's/+.*//'`
out_name=${tmp}_ANTS
print_out=${out_path}/${out_file}



if [[ ${#ANTSPATH} -le 3 ]] ; then
    setPath >&2
fi


if [ ! -s ${ANTSPATH}/antsRegistrationSyN.sh ]; then
    binCheck >&2
fi


argnum=$#

if [ $argnum == 0 ] || [ $argnum <= 2 ]; then
    Usage >&2
    exit 1
fi



echo "Starting antifyFuntional.sh for ${func_data}" > $print_out




#################################
#####  Run Checks            ####
#################################

echo "Performing checks:" >> $print_out

#Check that the deform matrix exists (and decompress it)
echo "Checking for ${deform_matrix}1Warp.nii.gz" >> $print_out

if [ -f ${deform_matrix}1Warp.nii.gz ]; then
    gzip -d ${deform_matrix}1Warp.nii.gz
fi

if [ ! -f ${deform_matrix}1Warp.nii ]; then
    echo "${deform_matrix}1Warp.nii not found!" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
    else
    echo "found" >> $print_out
fi


#Check for Affine.mat
echo "Checking for ${deform_matrix}0GenericAffine.mat" >> $print_out

if [ ! -f ${deform_matrix}0GenericAffine.mat ]; then
    echo "${deform_matrix}0GenericAffine.mat not found!" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
    else
    echo "found" >> $print_out
fi


#Check that the model exists
echo "Checking that $model exists" >> $print_out

if [ ! -f $model ]; then
    echo "${model} not found!" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
    else
    echo "found" >> $print_out
fi


#Check that the func data exists
echo "Checking for $func_data " >> $print_out

if [ ! -f ${func_data}.HEAD ]; then
    echo "${func_data}.HEAD not found!" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
    else
    echo "found" >> $print_out
fi




#################################
#####  ANTIFY                ####
#################################

echo "" >> $print_out
echo "Starting ANTIFY......." >> $print_out


# Step 1: Get the maximum index of func data sub-briks (total sub-briks minus 1)
echo "Step 1: Get the maximum index of func data sub-briks (total sub-briks minus 1)" >> $print_out
ndatasets=`3dinfo -nvi $func_data`
echo "Success" >> $print_out


# Step 2: de-oblique the functional dataset
echo "Step 2: de-oblique the functional dataset" >> $print_out
3drefit -deoblique $func_data
echo "Success" >> $print_out


# Step 3: make a temporary split of the functional data and have it in NIFTI format
echo "Step 3: make a temporary split of the functional data and have it in NIFTI format" >> $print_out

for i in $(seq 0 $ndatasets); do
    3dAFNItoNIFTI -prefix tmp_${i} ${func_data}[${i}]

    if [ -f tmp_${i}.nii ]; then
        echo "Success for $i" >> $print_out
        else
        echo "Something went wrong with $i" >> $print_out
        echo "Exiting ...." >> $print_out
        exit 1
    fi

done


# Step 4: Run WarpImageMultiTransform
echo "Step 4: Run WarpImageMultiTransform" >> $print_out

for i in $(seq 0 $ndatasets); do
    WarpImageMultiTransform 3 tmp_$i.nii tmpal_$i.nii ${deform_matrix}1Warp.nii ${deform_matrix}0GenericAffine.mat -R ${model}

    if [ -f tmpal_$i.nii ]; then
        echo "Success for $i" >> $print_out
        else
        echo "Something went wrong with ANTs alignment on $i" >> $print_out
        echo "Exiting ...." >> $print_out
        exit 1
    fi

done


# Step 5: Put it back together
echo "Step 5: Put it back together" >> $print_out
cmd="3dbucket -prefix ${out_name} "
for i in $(seq 0 $ndatasets); do
    cmd="$cmd -fbuc tmpal_${i}.nii "
done
`$cmd`

if [ -f ${out_name}+orig.HEAD ]; then
    echo "Success" >> $print_out
    else
    echo "Something went wrong with 3dbucket" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
fi


# Step 6: resample from 1x1x1mm to whatever the resolution used to be
echo "Step 6: resample from 1x1x1mm to whatever the resolution used to be" >> $print_out
#resolution=`3dAttribute DELTA $func_data`
resolution=`3dinfo -ad3 $func_data`
3dresample -dxyz $resolution -prefix ${out_name}_resampled+orig -inset ${out_name}+orig -rmode Cubic

if [ -f ${out_name}_resampled+orig.HEAD ]; then
    echo "Success" >> $print_out
    else
    echo "Something went wrong with resampling" >> $print_out
    echo "Exiting ...." >> $print_out
    exit 1
fi


# Step 7: Put the sub-brik labels back on
echo "Step 7: Put the sub-brik labels back on" >> $print_out

for i in $(seq 0 $ndatasets); do
    label=`3dinfo -label ${func_data}[${i}]`
    3drefit -sublabel $i ${label} ${out_name}_resampled+orig
done

echo "Success" >> $print_out


# Step 8: make it tlrc and MNI
echo "Step 8: make it tlrc and MNI" >> $print_out
3drefit -view tlrc -space MNI ${out_name}_resampled+orig

echo "Success" >> $print_out


# Step 9: clean up
rm tmp*.nii ${out_name}+orig*


# Step 10: Brag about it
cat >> $print_out << EOF
Finished antifyFunctional.sh for $func_data

##############################



EOF
