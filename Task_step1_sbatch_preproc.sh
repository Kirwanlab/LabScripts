#!/bin/bash

#SBATCH --time=30:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "TS1"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE




# Written by Nathan Muncy on 10/24/18

# --- Notes, in no particular order
#
# 1) will do first preproc steps - volreg, align, tshift, scale, mask construction
#
# 2) places needing researcher input are marked by "###??? update ..."
#
# 3) big steps are marked "### --- foobar --- ###", annotations with "#"
#
# 4) assumes dcm2niix/BIDS format
#
# 5) can accept blip scans (P->A), and should have one blip scan per phase




subj=$1


###??? update these
parDir=~/compute/STT_reml
workDir=${parDir}/derivatives/$subj
dataDir=${parDir}/rawdata/$subj

tempDir=~/bin/Templates/vold2_mni
template=${tempDir}/vold2_mni_brain+tlrc
priorDir=${tempDir}/priors_ACT

blip=0										# blip toggle (1=on)

phaseArr=(STUDY TEST{1,2})					# Each PHASE of experiment, in order (TEST1 precedes TEST2 but followed STUDY). This is absolutely necessary - so put something here
blockArr=(1 2 4)   							# number of blocks (runs) in each phase. E.g. STUDY had 1 block, TEST1 had 2 blocks. INTEGER!
phaseLen=${#phaseArr[@]}




### --- Set up --- ###
#
# Copies data to derivative $workDir, determines number
# of EPI runs ($block, $numRuns).
# Blocks are named according to their phase


### Checks
# check set up
cd ${dataDir}/func

if [ ${#phaseArr[@]} != ${#blockArr[@]} ]; then
	echo "\$phaseArr and \$blockArr are not same length. Exit 1" >&2
	exit 1
fi


# check blocks, names
for i in ${blockArr[@]}; do
	let totBlock+=$i
done

runCount=`ls *run*.nii.gz | wc -l`
if [ $runCount != $totBlock ]; then
	echo "Sum of \$blockArr and number of runs are not equal. Exit 2" >&2
	exit 2
fi

# check blips
if [ $blip == 1 ]; then
	blipCount=`ls *phase*.nii.gz | wc -l`
	if [ $blipCount != $phaseLen ]; then
		echo "Number of blips != number of phases. Exit 3" >&2
		exit 3
	fi
fi



### Copy data
# 3dcopy/rename epi data according to phase membership, determine number of blocks and set block arr
c=0; for i in *run*.nii.gz; do

	tmp=${i%_bold*}
	run=${tmp##*_}

	unset hold phaseName
	x=0; while [ $x -lt $phaseLen ] && [ -z $phaseName ]; do
		let hold+=${blockArr[$x]}
		if [ $(($c+1)) -le $hold ]; then
			phaseName=${phaseArr[$x]}
		fi
		let x=$[x+1]
	done

	if [ ! -f ${workDir}/${run}_${phaseName}+orig.HEAD ]; then
		3dcopy $i ${workDir}/${run}_${phaseName}+orig
	fi

	block[$c]=${run}_${phaseName}
	let c=$[$c+1]
done

numRuns=${#block[@]}


# blip data
if [ $blip == 1 ]; then
	for((i=1; i<=$blipCount; i++)); do
		if [ ! -f ${workDir}/phase${i}_Blip+orig.HEAD ]; then
			3dcopy *phase-${i}_blip.nii.gz ${workDir}/phase${i}_Blip+orig
		fi
	done
fi


# get t1 data
cd ${dataDir}/anat

if  [ ! -f ${workDir}/struct+orig.HEAD ]; then
	3dcopy *_T1w.nii.gz ${workDir}/struct+orig
fi




### --- Blip --- ###
#
# Correct for signal fallout in OFC by using data from an opposite phase encoding (Rev).
# Find median datasets, compute midpoints, warp median datasets,
# warp EPI timeseries


cd $workDir

if [ $blip == 1 ]; then
	if [ ! -f run-1_${phaseArr[0]}_blip+orig.HEAD ]; then
		for a in ${!phaseArr[@]}; do

			# median blip
			num=$(($a+1))
			phase=${phaseArr[$a]}

			3dTcat -prefix tmp_blip${num} phase${num}_Blip+orig
			3dTstat -median -prefix tmp_blip${num}_med tmp_blip${num}+orig
			3dAutomask -apply_prefix tmp_blip${num}_med_masked tmp_blip${num}_med+orig


			# median task
			unset tcat_input
			for b in run-*${phase}+orig.HEAD; do
				tcat_input+="${b%.*} "
			done

			3dTcat -prefix tmp_all${phase} $tcat_input
			3dTstat -median -prefix tmp_all${phase}_med tmp_all${phase}+orig
			3dAutomask -apply_prefix tmp_all${phase}_med_masked tmp_all${phase}_med+orig


			# compute midpoint
			3dQwarp -plusminus -pmNAMES Rev Task \
				-pblur 0.05 0.05 -blur -1 -1 \
				-noweight -minpatch 9 \
				-source tmp_blip${num}_med_masked+orig \
				-base tmp_all${phase}_med_masked+orig \
				-prefix tmp_${phase}_warp


			# warp median dsets/masks
			3dNwarpApply -quintic -nwarp tmp_${phase}_warp_Task_WARP+orig \
				-source tmp_all${phase}_med+orig \
				-prefix tmp_all${phase}_NWA

			3drefit -atrcopy tmp_all${phase}+orig IJK_TO_DICOM_REAL tmp_all${phase}_NWA+orig

			3dNwarpApply -quintic -nwarp tmp_${phase}_warp_Task_WARP+orig \
				-source tmp_all${phase}_med_masked+orig \
				-prefix tmp_all${phase}_med_Task_masked

			3drefit -atrcopy tmp_all${phase}+orig IJK_TO_DICOM_REAL tmp_all${phase}_med_Task_masked+orig


			# warp EPI
			for j in run-*_${phase}+orig.HEAD; do

				scan=${j%+*}
				if [ ! -f ${scan}_blip+orig.HEAD ]; then

					3dNwarpApply -quintic -nwarp tmp_${phase}_warp_Task_WARP+orig \
						-source ${j%.*} \
						-prefix ${scan}_blip

					3drefit -atrcopy tmp_all${phase}+orig IJK_TO_DICOM_REAL ${scan}_blip+orig
				fi
			done
			rm tmp_*
		done
	fi
fi





### --- Volreg Setup --- ###
#
# Outliers will be detected for later exclusion. The TR of the
# experiment with the minimum noise will become the volume registration
# base and used to construct a volreg_base file.
#
# The blip file will now be incorporated to account for signal fallout


for j in ${block[@]}; do

	if [ $blip == 1 ]; then
		input=${j}_blip+orig
	else
		input=${j}+orig
	fi

	# build outcount list
	hold=`3dinfo -ntimes $input`
	tr_counts+="$hold "


	if [ ! -s outcount.${j}.1D ]; then

		# determine polort arg
		len_tr=`3dinfo -tr $input`
		pol_time=$(echo $(echo $hold*$len_tr | bc)/150 | bc -l)
		pol=$((1 + `printf "%.0f" $pol_time`))

		3dToutcount -automask -fraction -polort $pol -legendre $input > outcount.${j}.1D


		# censor
		> out.${j}.pre_ss_warn.txt
		1deval -a outcount.${j}.1D -expr "1-step(a-0.1)" > out.cen.${j}.1D
		if [ `1deval -a outcount.${j}.1D"{0}" -expr "step(a-0.4)"` ]; then
			echo "** TR #0 outliers: possible pre-steady state TRs in run $j"  >> out.${j}.pre_ss_warn.txt
		fi
	fi
done


if [ ! -f epi_vr_base+orig.HEAD ]; then

	# determine min volume
	cat outcount.*.1D > outcount_all.1D
	minindex=`3dTstat -argmin -prefix - outcount_all.1D\'`
	ovals=(`1d_tool.py -set_run_lengths $tr_counts -index_to_run_tr $minindex`)

	minoutrun=${ovals[0]}
	minouttr=${ovals[1]}


	### Time shift here, skipping for multiband/dcm2nii.
	# T-shifted base would be input for next step (3dbucket)


	# determine volreg base by matching $ovals, could be in any run
	c=0; for ((d=1; d <= $numRuns; d++)); do
		if [ 0$d == $minoutrun ]; then
			baseRun=${block[$c]}
		fi
		let c=$[$c+1]
	done


	# construct volreg base, print out
	if [ $blip == 1 ]; then
		newBase=${baseRun}_blip
	else
		newBase=$baseRun
	fi

	3dbucket -prefix epi_vr_base ${newBase}+orig"[${minouttr}]"
	echo "$minoutrun $minouttr $newBase" > out_vr_base.txt
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


if [ ! -s anat.un.aff.Xat.1D ]; then

	# calc align of epi/anat
	align_epi_anat.py \
	-anat2epi \
	-anat struct+orig \
	-save_skullstrip \
	-suffix _al_junk \
	-epi epi_vr_base+orig \
	-epi_base 0 \
	-epi_strip 3dAutomask \
	-cost lpc+ZZ \
	-volreg off \
	-tshift off


	# calc non-linear warp
	auto_warp.py -base $template -input struct_ns+orig -skull_strip_input no
	3dbucket -prefix struct_ns awpy/struct_ns.aw.nii*
	cp awpy/anat.un.aff.Xat.1D .
	cp awpy/anat.un.aff.qw_WARP.nii .
fi


# determine voxel size - assumes isotropic voxels
gridSize=`3dinfo -di ${block[0]}+orig`


for j in ${block[@]}; do
	if [ ! -f tmp_${j}_mask_warped+tlrc.HEAD ]; then

		if [ $blip == 1 ]; then
			input=${j}_blip+orig
		else
			input=${j}+orig
		fi

		# calc volreg
		3dvolreg -verbose -zpad 1 -base epi_vr_base+orig \
		-1Dfile dfile.${j}.1D -prefix ${j}_volreg \
		-cubic \
		-1Dmatrix_save mat.${j}.vr.aff12.1D \
		$input


		# concat calcs for epi movement (volreg, align, warp)
		cat_matvec -ONELINE \
		anat.un.aff.Xat.1D \
		struct_al_junk_mat.aff12.1D -I \
		mat.${j}.vr.aff12.1D > mat.${j}.warp.aff12.1D


		# warp epi
		3dNwarpApply -master struct_ns+tlrc \
		-dxyz $gridSize \
		-source $input \
		-nwarp "anat.un.aff.qw_WARP.nii mat.${j}.warp.aff12.1D" \
		-prefix tmp_${j}_nomask


		# warp mask for extents masking; make intersection mask (epi+anat)
		3dcalc -overwrite -a $input -expr 1 -prefix tmp_${j}_mask

		3dNwarpApply -master struct_ns+tlrc \
		-dxyz $gridSize \
		-source tmp_${j}_mask+orig \
		-nwarp "anat.un.aff.qw_WARP.nii mat.${j}.warp.aff12.1D" \
		-interp cubic \
		-ainterp NN -quiet \
		-prefix tmp_${j}_mask_warped

		3dTstat -min -prefix tmp_${j}_min tmp_${j}_mask_warped+tlrc
	fi
done


# create extents mask, delete TRs w/missing data
for j in ${block[@]}; do
	if [ ! -f ${j}_volreg_clean+tlrc.HEAD ]; then

		if [ $numRuns > 1 ]; then
			3dMean -datum short -prefix tmp_mean_${j} tmp_${j}*_min+tlrc.HEAD
			3dcalc -a tmp_mean_${j}+tlrc -expr 'step(a-0.999)' -prefix ${j}_epiExt_mask
		else
			3dcopy tmp_${j}_min+tlrc.HEAD ${j}_epiExt_mask
		fi

		3dcalc -a tmp_${j}_nomask+tlrc -b ${j}_epiExt_mask+tlrc -expr 'a*b' -prefix ${j}_volreg_clean
	fi
done


# warp volreg base
if [ ! -f final_epi_vr_base+tlrc.HEAD ]; then

	# concat align, warp calcs
	cat_matvec -ONELINE \
	anat.un.aff.Xat.1D \
	struct_al_junk_mat.aff12.1D -I  > mat.basewarp.aff12.1D

	3dNwarpApply -master struct_ns+tlrc \
	-dxyz $gridSize \
	-source epi_vr_base+orig \
	-nwarp "anat.un.aff.qw_WARP.nii mat.basewarp.aff12.1D" \
	-prefix final_epi_vr_base
fi


# anat copy
if [ ! -f final_anat+tlrc.HEAD ]; then
	3dcopy struct_ns+tlrc final_anat
fi


# record registration costs; affine warp follower dsets
if [ ! -f final_anat_head+tlrc.HEAD ]; then

	3dAllineate -base final_epi_vr_base+tlrc -allcostX  \
	-input final_anat+tlrc | tee out.allcostX.txt

	3dNwarpApply -source struct+orig \
	-master final_anat+tlrc \
	-ainterp wsinc5 \
	-nwarp anat.un.aff.qw_WARP.nii anat.un.aff.Xat.1D \
	-prefix final_anat_head
fi




### --- Create Masks --- ###
#
# An EPI T1 intersection mask is constructed, then tissue-class
# masks are created (these are used for REML). The AFNI
# version of tiss-seg is left, but I prefer the Atropos priors.


# union inputs (combine Run masks); anat mask; intersecting; group
if [ ! -f final_anat_mask+tlrc.HEAD ]; then

	for j in ${block[@]}; do
		3dAutomask -prefix tmp_mask.${j} ${j}_volreg_clean+tlrc
	done
	3dmask_tool -inputs tmp_mask.*+tlrc.HEAD -union -prefix full_mask

	3dresample -master full_mask+tlrc -input struct_ns+tlrc -prefix tmp_anat_resamp
	3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_anat_resamp+tlrc -prefix final_anat_mask

	3dmask_tool -input full_mask+tlrc final_anat_mask+tlrc -inter -prefix mask_epi_anat
	3dABoverlap -no_automask full_mask+tlrc final_anat_mask+tlrc | tee out.mask_ae_overlap.txt

	3dresample -master full_mask+tlrc -prefix ./tmp_resam_group -input $template
	3dmask_tool -dilate_input 5 -5 -fill_holes -input tmp_resam_group+tlrc -prefix Template_mask
fi


### segment tissue class, generate masks

## AFNI way
#3dSeg -anat final_anat+tlrc -mask AUTO -classes 'CSF ; GM ; WM'

#for j in CSF GM WM; do

	#3dmask_tool -input Segsy/Classes+tlrc"<${j}>" -prefix tmp_mask_$j
	#3dresample -master ${block[0]}_volreg_clean+tlrc -rmode NN -input tmp_mask_${j}+tlrc -prefix final_mask_${j}
	#3dmask_tool -input Segsy/Classes+tlrc"<${j}>" -dilate_input -1 -prefix tmp_mask_${j}_eroded
	#3dresample -master ${block[0]}_volreg_clean+tlrc -rmode NN -input tmp_mask_${j}_eroded+tlrc -prefix final_mask_${j}_eroded
#done


# seg tissue class, with Atropos priors, for REML step
if [ ! -f final_mask_GM_eroded+tlrc.HEAD ]; then

	# get priors
	tiss=(CSF GMc WM GMs)
	prior=(Prior{1..4})
	tissN=${#tiss[@]}

	c=0; while [ $c -lt $tissN ]; do
		cp ${priorDir}/${prior[$c]}.nii.gz ./tmp_${tiss[$c]}.nii.gz
		let c=$[$c+1]
	done
	c3d tmp_GMc.nii.gz tmp_GMs.nii.gz -add -o tmp_GM.nii.gz

	# resample, erode
	for i in CSF GM WM; do

		c3d tmp_${i}.nii.gz -thresh 0.3 1 1 0 -o tmp_${i}_bin.nii.gz
		3dresample -master ${block[0]}_volreg_clean+tlrc -rmode NN -input tmp_${i}_bin.nii.gz -prefix final_mask_${i}+tlrc
		3dmask_tool -input tmp_${i}_bin.nii.gz -dilate_input -1 -prefix tmp_mask_${i}_eroded
		3dresample -master ${block[0]}_volreg_clean+tlrc -rmode NN -input tmp_mask_${i}_eroded+orig -prefix final_mask_${i}_eroded
	done
fi





### --- Scale --- ###
#
# Data is scaled by mean signal - gotta reduce them confounds.


for j in ${block[@]}; do
	if [ ! -f ${j}_scale+tlrc.HEAD ]; then

		3dTstat -prefix tmp_tstat_$j ${j}_volreg_clean+tlrc

		3dcalc \
		-a ${j}_volreg_clean+tlrc \
		-b tmp_tstat_${j}+tlrc \
		-c ${j}_epiExt_mask+tlrc \
		-expr 'c * min(200, a/b*100)*step(a)*step(b)' \
		-prefix ${j}_scale
	fi
done
