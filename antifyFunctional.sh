#!/bin/bash



function Usage {
cat << USAGE

##############################

 antifyFunctional.sh

 Script to apply ANTS transformation to afni functional data.

 Original perl script by Craig Stark
 Modified by Brock Kirwan

 Update 10/12/2015: Added affine transformation to the WarpImageMultiTransform call
                    Added "deoblique" step before the afni to nii conversion to deal with alignment issues
        11/19/2016: Ported to bash

 -----------------------------

 Usage:  antifyFunctional.sh deformation_matrix model functional_dataset

 antifyFunctional warps functional data into template space, using calculations previously performed
 by registering T1-structural files to a template. It requires both the Warp.nii and Affine.txt
 from ANTs registration.

 Required Arguments:

    'deformation_matrix' is the prefix of the deformation matrix file (usually struct_rotated).
    'model' is the fullpath to your model template.
    'functional_dataset' is the prefix of the functional dataset.

##############################

USAGE
    exit 1
}




NUMPARAMS=$#

if [ ${#ANTSPATH} -le 3 ] ; then
echo we guess at your ants path
export ANTSPATH=${ANTSPATH:="/usr/local/bin/ants"} # a guess based on my machine
fi
if [ ! -s ${ANTSPATH}/ANTS ]; then
echo we cant find the ANTS program -- does not seem to exist.  please \(re\)define \$ANTSPATH in your environment.
exit
fi


if [ $NUMPARAMS -lt 3  ]; then
    Usage >&2
    exit 1
fi

#Parse the command line input
deform_matrix=$1
model=$2
func_data=$3
tmp=`echo $func_data | sed 's/+.*//'`
out_name=${tmp}_ANTS

#Check that the deform matrix exists (and decompress it)
if [ -f ${deform_matrix}Warp.nii.gz ]; then
gzip -d ${deform_matrix}Warp.nii.gz
fi

if [ ! -f ${deform_matrix}Warp.nii ]; then
echo "${deform_matrix}Warp.nii not found!"
exit
fi

#Check that the model exists
if [ ! -f $model ]; then
echo "${model} not found!"
exit
fi

#Check that the func data exists
if [ ! -f ${func_data}.HEAD ]; then
echo "${func_data}.HEAD not found!"
exit
fi

# Step 1: Get the maximum index of func data sub-briks (total sub-briks minus 1)
ndatasets=`3dinfo -nvi $func_data`

# Step 2: de-oblique the functional dataset
echo "Making the functional dataset de-oblique"
3dRefit -deoblique $func_data

# Step 3: make a temporary split of the functional data and have it in NIFTI format
echo
for i in $(seq 0 $ndatasets); do
echo "."
3dAFNItoNIFTI -prefix tmp_${i} ${func_data}[${i}]
done

# Step 4: Run WarpImageMultiTransform
echo "Running alighment"
for i in $(seq 0 $ndatasets); do
WarpImageMultiTransform 3 tmp_$i.nii tmpal_$i.nii ${deform_matrix}Warp.nii -R ${model} ${deform_matrix}Affine.txt
if [ ! -f tmpal_$i.nii ]; then
echo "Something went wrong with ANTs alignment"
exit
fi
done

# Step 5: Put it back together
echo "Reassembling"
cmd="3dBucket -prefix ${out_name} "
for i in $(seq 0 $ndatasets); do
cmd="$cmd -fbuc tmpal_${i}.nii "
done
`$cmd`

# Step 6: resample from 1x1x1mm to whatever the resolution used to be
echo "Resampling"
resolution=`3dAttribute DELTA ${func_data}`
#force it to be positive
resolution=`echo $resolution | sed 's/-//g'`
3dresample -dxyz $resolution -prefix ${out_name}_resampled+orig -inset ${out_name}+orig -rmode Cubic

# Step 7: Put the sub-brik labels back on
echo "Relabling"
for i in $(seq 0 $ndatasets); do
label=`3dinfo -label ${func_data}[${i}]`
3drefit -sublabel $i ${label} ${out_name}_resampled+orig
done

# Step 8: make it tlrc and MNI
3drefit -view tlrc -space MNI ${out_name}_resampled+orig

# Step 9: clean up
rm tmp*.nii ${out_name}+orig*

# Step 10: Brag about it
echo "Finished antifyFunctional.sh"
