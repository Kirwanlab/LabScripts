#!/bin/bash

#SBATCH --time=20:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "RS1"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE




# written by Nathan Muncy on 12/3/19


###--- Notes
#
# This script does basic pre-processing for resting state data
#
# Major steps following the comment ### --- Foo --- ###,
# minor steps are also annotated
#
# Sections ###??? require research attention
#
# Phase = session
#
# Researcher Input A) location of parent directory
# Researcher Input B) a bold file in rawdata
# Researcher Input C) verify name of bold file in rawdata
# Researcher Input D) this was written for a sleep study, where a T1-weighted file only exists for one session. Update accordingly





subj=$1
phase=$2	# required


### --- Set Variables --- ###

parDir=~/compute/MNI_NREM/MRI_processed								###??? Research Input A
workDir=${parDir}/derivatives/${subj}/ses-$phase 					# subject, session directory
rawDir=${parDir}/rawdata/$subj 										# location of subject raw data
refFile=${parDir}/rawdata/sub-021/ses-Sleep/func/sub-021_RS_run-Sleep_bold.nii.gz		###??? Research Input B

blip=0																# if blip data exists (1=on)
multiBand=1															# if vibration is present from MB protocol (1=on)

# Location of template/priors
tempDir=~/bin/Templates/vold2_mni
template=${tempDir}/vold2_mni_brain+tlrc
actDir=${tempDir}/priors_ACT




### --- Set up --- ###
#
# Run checks to make sure the script is set up correctly,
# and copy data from rawdata


# Check
for i in {1..4}; do
	if [ ! -f ${actDir}/Prior${i}.nii.gz ]; then
		echo >&2
		echo "Prior${i} not detected. Check \$actDir. Exit 1" >&2
		echo >&2
		exit 1
	fi
done

if [ ! -f ${template}.HEAD ]; then
	echo >&2
	echo "Template not detected. Check \$tempDir. Exit 2" >&2
	echo >&2
	exit 2
fi


# Copy Data
mkdir -p $workDir

if [ ! -f ${workDir}/Resting_${phase}+orig.HEAD ]; then
	3dcopy ${rawDir}/ses-${phase}/func/${subj}_RS_run-${phase}_bold.nii.gz ${workDir}/Resting_${phase}+orig 	###??? Researcher Input C
fi

if [ ! -f ${workDir}/struct+orig.HEAD ]; then
	3dcopy ${rawDir}/ses-Sleep/anat/${subj}_T1w.nii.gz ${workDir}/struct+orig 					###??? Researcher Input D

fi

if [ ! -f ${workDir}/Resting_${phase}+orig.HEAD ] || [ ! -f ${workDir}/struct+orig.HEAD ]; then
	echo >&2
	echo "Func and/or Anat data not detected in $workDir. Check $rawDir. Exit 3" >&2
	echo >&2
	exit 3
fi




### --- Blip --- ###
#
# Correct for signal fallout in OFC by using data from an opposite phase encoding (Rev).
# Find median datasets, compute midpoints, warp median datasets,
# warp EPI timeseries


cd $workDir

# build outcount list
tr_counts=`3dinfo -ntimes Resting_${phase}+orig`
gridSize=`3dinfo -di Resting_${phase}+orig`
checkSize=`3dinfo -di $refFile`


## patch for FOV problem: verify voxel dim, resample if needed
if [ $gridSize != $checkSize ]; then

	3dcopy Resting_${phase}+orig Resting_origRes_${phase} && rm Resting_${phase}+orig*
	3dresample -master $refFile -rmode NN -input Resting_origRes_${phase}+orig -prefix Resting_${phase}
	gridSize=`3dinfo -di Resting_${phase}+orig`
fi

# check
if [ $gridSize != $checkSize ]; then
	echo >&2
	echo "Voxel dimensionality not correct. Exiting 4" >&2
	echo >&2; exit 4
fi


# do blip
if [ $blip == 1 ]; then
	if [ ! -f Resting_${phase}_noBlip+orig.HEAD ]; then

		# medain blip
		3dTcat -prefix tmp_blip ${rawDir}/ses-${phase}/func/${subj}_task-restblip_bold.nii.gz
		3dTstat -median -prefix tmp_blip_med tmp_blip+orig
		3dAutomask -apply_prefix tmp_blip_med_masked tmp_blip_med+orig


		# median rs epi
		3dcopy Resting_${phase}+orig Resting_${phase}_noBlip+orig && rm Resting_${phase}+orig.*

		3dTcat -prefix tmp_all Resting_${phase}_noBlip+orig
		3dTstat -median -prefix tmp_all_med tmp_all+orig
		3dAutomask -apply_prefix tmp_all_med_masked tmp_all_med+orig


		# compute midpoint
		3dQwarp -plusminus -pmNAMES Rev Task \
			-pblur 0.05 0.05 -blur -1 -1 \
			-noweight -minpatch 9 \
			-source tmp_blip_med_masked+orig \
			-base tmp_all_med_masked+orig \
			-prefix tmp_warp


		# warp median dsets/masks
		3dNwarpApply -quintic -nwarp tmp_warp_Task_WARP+orig \
			-source tmp_all_med+orig \
			-prefix tmp_all_NWA

		3drefit -atrcopy tmp_all+orig IJK_TO_DICOM_REAL tmp_all_NWA+orig

		3dNwarpApply -quintic -nwarp tmp_warp_Task_WARP+orig \
			-source tmp_all_med_masked+orig \
			-prefix tmp_all_med_Task_masked

		3drefit -atrcopy tmp_all+orig IJK_TO_DICOM_REAL tmp_all_med_Task_masked+orig


		# warp EPI
		3dNwarpApply -quintic -nwarp tmp_warp_Task_WARP+orig \
			-source Resting_${phase}_noBlip+orig \
			-prefix Resting_${phase}+orig

		3drefit -atrcopy tmp_all+orig IJK_TO_DICOM_REAL Resting_${phase}+orig
	fi
fi





### --- Volreg Setup --- ###
#
# Outliers will be detected for later exclusion. The volume of the
# experiment with the minimum noise will be extracted and serve
# as volume registration base.


if [ ! -s outcount.RS.1D ]; then

	# determine polort arg, outliers
	len_tr=`3dinfo -tr Resting_${phase}+orig`
	pol_time=$(echo $(echo $tr_counts*$len_tr | bc)/150 | bc -l)
	pol=$((1 + `printf "%.0f" $pol_time`))
	3dToutcount -automask -fraction -polort $pol -legendre Resting_${phase}+orig > outcount.RS.1D

	# censor - more conservative threshold (0.05) used for RS
	> out.RS.pre_ss_warn.txt
	1deval -a outcount.RS.1D -expr "1-step(a-0.05)" > out.cen.RS.1D
	if [ `1deval -a outcount.RS.1D"{0}" -expr "step(a-0.4)"` ]; then
		echo "** TR #0 outliers: possible pre-steady state TRs in run RS"  >> out.RS.pre_ss_warn.txt
	fi
fi


# Despike, get min outlier volume
if [ ! -f vr_base+orig.HEAD ]; then

    3dDespike -NEW -nomask -prefix tmp_despike Resting_${phase}+orig

	minindex=`3dTstat -argmin -prefix - outcount.RS.1D\'`
	ovals=(`1d_tool.py -set_run_lengths $tr_counts -index_to_run_tr $minindex`)
	minouttr=${ovals[1]}
    3dbucket -prefix vr_base tmp_despike+orig"[$minouttr]"
fi




### --- Normalize Data --- ###
#
# First, a rigid transformation with a function will be calculated
# bx epi & t1. Skull-stripping happens in this step. Second a
# non-linear diffeomorphich transformation of rotated brain to
# template space is calculated. Third, we get the volreg calculation.
# EPI is warped into template space with a single interpolation, by
# combining the rigid, volreg, and diffeo calculations. T1 also warped,
# as is the volreg_base by using the appropriate calcs. An extents
# mask is constructed and used to delete TRs with missing data.
# Registration cost is recorded


# align
if [ ! -f struct_ns+orig.HEAD ]; then

	align_epi_anat.py \
	-anat2epi \
	-anat struct+orig \
	-save_skullstrip \
	-suffix _rotated \
	-epi vr_base+orig \
	-epi_base 0 \
	-epi_strip 3dAutomask \
	-cost lpc+ZZ \
	-volreg off \
	-tshift off
fi


# normalize
if [ ! -f anat.un.aff.qw_WARP.nii ]; then
	auto_warp.py -base $template -input struct_ns+orig -skull_strip_input no
	3dbucket -prefix struct_ns awpy/struct_ns.aw.nii*
	cp awpy/anat.un.aff.Xat.1D .
	cp awpy/anat.un.aff.qw_WARP.nii .
fi


if [ ! -f struct_ns+tlrc.HEAD ]; then
	echo >&2
	echo "Normalization failed - no struct_ns+tlrc.HEAD detected. Exit 5" >&2
	echo >&2; exit 5
fi


### Move data to template space
#
# 3 calcs are used - 1) volreg, 2) rotating struct,
# and 3) normalization. Epi data is moved from native
# into template space via 1+(2^-1)+3.
#
# Patch for multiband - using volreg to account for
# vibration in sequence. Move volreg data to template
# via (2^-1)+3.

if [ ! -f tmp_epi_mask_warped+tlrc.HEAD ]; then

	# calc volreg
	3dvolreg -verbose -zpad 1 -base vr_base+orig \
	-1Dfile dfile_Resting_${phase}.1D -prefix tmp_epi_volreg \
	-cubic \
	-1Dmatrix_save mat.vr.aff12.1D \
	tmp_despike+orig


	# concat calcs for epi movement (volreg, align, warp)
	if [ $multiBand == 1 ]; then
		cat_matvec -ONELINE \
		anat.un.aff.Xat.1D \
		struct_rotated_mat.aff12.1D -I > mat.warp_epi.aff12.1D
	else
		cat_matvec -ONELINE \
		anat.un.aff.Xat.1D \
		struct_rotated_mat.aff12.1D -I \
		mat.vr.aff12.1D > mat.warp_epi.aff12.1D
	fi


	# warp epi
	if [ $multiBand == 1 ]; then
		epiMov=tmp_epi_volreg+orig
	else
		epiMov=tmp_despike+orig
	fi

	3dNwarpApply -master struct_ns+tlrc \
	-dxyz $gridSize \
	-source $epiMov \
	-nwarp "anat.un.aff.qw_WARP.nii mat.warp_epi.aff12.1D" \
	-prefix tmp_epi_warped_nomask


	# warp mask for extents masking; make intersection mask (epi+anat)
	3dcalc -overwrite -a $epiMov -expr 1 -prefix tmp_epi_mask

	3dNwarpApply -master struct_ns+tlrc \
	-dxyz $gridSize \
	-source tmp_epi_mask+orig \
	-nwarp "anat.un.aff.qw_WARP.nii mat.warp_epi.aff12.1D" \
	-interp cubic \
	-ainterp NN -quiet \
	-prefix tmp_epi_mask_warped
fi


# create epi extents mask
if [ ! -f tmp_epi_clean+tlrc.HEAD ]; then

	3dTstat -min -prefix tmp_epi_min tmp_epi_mask_warped+tlrc
	3dcopy tmp_epi_min+tlrc mask_epi_extents
	3dcalc -a tmp_epi_warped_nomask+tlrc -b mask_epi_extents+tlrc -expr 'a*b' -prefix tmp_epi_clean
fi


# warp volreg base into template space
if [ ! -f final_vr_base+tlrc.HEAD ];then

	# concat align, warp calcs
	cat_matvec -ONELINE \
	anat.un.aff.Xat.1D \
	struct_rotated_mat.aff12.1D -I  > mat.basewarp.aff12.1D

	3dNwarpApply -master struct_ns+tlrc \
	-dxyz $gridSize \
	-source vr_base+orig \
	-nwarp "anat.un.aff.qw_WARP.nii mat.basewarp.aff12.1D" \
	-prefix final_vr_base
fi


# anat copy
if [ ! -f final_anat+tlrc.HEAD ]; then
	3dcopy struct_ns+tlrc final_anat
fi


# record registration costs; affine warp follower dsets
if [ ! -f final_anat_head+tlrc.HEAD ]; then

	3dAllineate -base final_vr_base+tlrc -allcostX -input final_anat+tlrc | tee out.allcostX.txt

	3dNwarpApply \
	-source struct+orig \
	-master final_anat+tlrc \
	-ainterp wsinc5 \
	-nwarp anat.un.aff.qw_WARP.nii anat.un.aff.Xat.1D \
	-prefix final_anat_head
fi


# Blur data
if [ ! -f tmp_epi_blur+tlrc.HEAD ]; then

	int=`printf "%.0f" $gridSize`
	blur=$((2*$int))
	3dmerge -1blur_fwhm $blur -doall -prefix tmp_epi_blur tmp_epi_clean+tlrc
fi




### --- Create Masks --- ###
#
# An EPI T1 intersection mask is constructed, then tissue-class
# masks are created (these are used for REML).
# Tissue seg is based on Atropos Priors


# union inputs (combine Run masks); anat mask; intersecting; group
if [ ! -f final_anat_mask+tlrc.HEAD ]; then

	3dAutomask -prefix tmp_mask_epi tmp_epi_blur+tlrc
	3dmask_tool -inputs tmp_mask_epi+tlrc.HEAD -union -prefix full_mask

	3dresample -master full_mask+tlrc -input struct_ns+tlrc -prefix tmp_anat_resamp
	3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_anat_resamp+tlrc -prefix final_anat_mask

	3dmask_tool -input full_mask+tlrc final_anat_mask+tlrc -inter -prefix mask_epi_anat
	3dABoverlap -no_automask full_mask+tlrc final_anat_mask+tlrc | tee out.mask_ae_overlap.txt

	3dresample -master full_mask+tlrc -prefix ./tmp_resam_group -input $template
	3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_resam_group+tlrc -prefix Template_mask
fi


if [ ! -f final_mask_GM_eroded+tlrc.HEAD ]; then

	# get priors
	tiss=(CSF GMc WM GMs)
	prior=(Prior{1..4})
	tissN=${#tiss[@]}

	c=0; while [ $c -lt $tissN ]; do
		cp ${actDir}/${prior[$c]}.nii.gz ./tmp_${tiss[$c]}.nii.gz
		let c=$[$c+1]
	done
	c3d tmp_GMc.nii.gz tmp_GMs.nii.gz -add -o tmp_GM.nii.gz

	# resample, erode
	for i in CSF GM WM; do
		c3d tmp_${i}.nii.gz -thresh 0.3 1 1 0 -o tmp_${i}_bin.nii.gz
		3dresample -master tmp_epi_blur+tlrc -rmode NN -input tmp_${i}_bin.nii.gz -prefix tmp_mask_${i}+tlrc
		3dmask_tool -input tmp_${i}_bin.nii.gz -dilate_input -1 -prefix tmp_mask_${i}_eroded
		3dresample -master tmp_epi_blur+tlrc -rmode NN -input tmp_mask_${i}_eroded+orig -prefix final_mask_${i}_eroded
	done
fi




### --- Scale Data --- ###
#
# Scale time series by mean signal.


if [ ! -f Resting_${phase}_scale+tlrc.HEAD ]; then

	3dTstat -prefix tmp_tstat tmp_epi_blur+tlrc

	3dcalc \
	-a tmp_epi_blur+tlrc \
	-b tmp_tstat+tlrc \
	-c mask_epi_extents+tlrc \
	-expr 'c * min(200, a/b*100)*step(a)*step(b)' \
	-prefix Resting_${phase}_scale
fi
