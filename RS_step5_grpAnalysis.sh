#!/bin/bash

#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=10   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10gb   # memory per CPU core
#SBATCH -J "RS5"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE





# written by Nathan Muncy on 12/3/19


### --- Notes
#
# This script will create a gray matter, intersection mask and
# 	then generate and run MVM scripts. A Monte Carlo
#	simulation will be conducted.
#
# Note - the MVM script is hardcoded for the Sleep study. This
# 	script will build the input for the dataTable argument. The
#	rest should be written.
#
# Researcher Input A) Array of sessions/phases
# Researcher Input B) Reference directory, for building lists of files
# Researcher Input C) Reference bold file, for resampling
# Researcher Input D) List with columns for subject number, group membership, etc.
# Researcher Input E) Array of subjects
# Researcher Input F) Array of group membership for e/participant (categorical - between-subject factor)
# Researcher Input G) Input string for each participant in dataTable. Written for Sleep study, so update accordingly to your study
# Researcher Input H) Update MVM arguments





### --- Set up --- ###

parDir=~/compute/MNI_NREM/MRI_processed						###??? experiment directory
workDir=${parDir}/derivatives
outDir=${parDir}/Analyses/grpAnalysis

phaseArr=(Sleep Wake)										###??? Researcher input A
refDir=${workDir}/sub-021/ses-Sleep							###??? Researcher input B
refFile=${refDir}/Resting_Sleep_scale+tlrc					###??? Researcher input C

tempDir=~/bin/Templates/vold2_mni
priorDir=${tempDir}/priors_ACT
mask=Intersection_GM_mask+tlrc

refList=${parDir}/docs/subj_list.txt						###??? Researcher input D
arrSubj=(`cat $refList | awk '{print $1}'`)					###??? Researcher input E
arrGroup=(`cat $refList | awk '{print $2}'`)				###??? Researcher input F




### --- Functions --- ###

# search array for string
MatchString () {
	local e match="$1"
	shift
	for e; do
		[[ "$e" == "$match" ]] && return 0
	done
	return 1
}




### --- Make Lists --- ###

# ZTrans files
cd $refDir

c=0; for i in ZTrans*HEAD; do
	tmp=${i%+*}
	zSeeds[$c]=${tmp#*_}
	let c=$[$c+1]
done

if [ -z ${zSeeds[@]} ]; then
	echo >&2
	echo "List of z-trans files empty. Do they exist in $refDir ?. Exit 1" >&2
	echo >&2; exit 1
fi


## subject list
#
# Just using Sleep lists for now (who to exclude - arrRem).
# Set arrays for subject and group, account for excluded participants

arrRem=(`cat ${outDir}/info_${phaseArr[0]}_rmSubj.txt`)
unset subjList groupList

c=0; for i in ${!arrSubj[@]}; do
	subj=sub-${arrSubj[$i]/P}
	if [ -d ${workDir}/${subj}/ses-${phaseArr[0]} ]; then
		MatchString "$subj" "${arrRem[@]}"
		if [ $? == 1 ]; then
			subjList[$c]=$subj
			groupList[$c]=${arrGroup[$i]}
			let c=$[$c+1]
		fi
	fi
done

# check
if [ ${#subjList[@]} != ${#groupList[@]} ]; then
	echo >&2
	echo "Subject and group lists are not same length. Exit 2" >&2
	echo >&2; exit 2
fi




### --- Create Masks --- ###
#
# This section will create a group mean intersection mask
# then threshold it at $thr to create a binary intersection mask.
# A gray matter mask will be constructed, and then the GM mask
# will be multiplied with the intersection mask to create a
# single GM intersection mask


# intersection mask
cd $outDir

if [ ! -f Group_epi_mask.nii.gz ]; then

	unset list
	for i in ${subjList[@]}; do
		for j in ${phaseArr[@]}; do
			list+="${workDir}/${i}/ses-${j}/mask_epi_anat+tlrc "
		done
	done

	3dMean -prefix ${outDir}/Group_epi_mean.nii.gz $list
	3dmask_tool -input $list -frac 0.3 -prefix ${outDir}/Group_epi_mask.nii.gz
fi


# make $mask
if [ ! -f ${mask}.HEAD ]; then

	# GM mask
	c3d ${priorDir}/Prior2.nii.gz ${priorDir}/Prior4.nii.gz -add -o tmp_Prior_GM.nii.gz
	3dresample -master $refFile -rmode NN -input tmp_Prior_GM.nii.gz -prefix tmp_Template_GM_mask.nii.gz

	# combine GM and intersection mask
	c3d tmp_Template_GM_mask.nii.gz Group_epi_mask.nii.gz -multiply -o tmp_Intersection_GM_prob_mask.nii.gz
	c3d tmp_Intersection_GM_prob_mask.nii.gz -thresh 0.1 1 1 0 -o tmp_Intersection_GM_mask.nii.gz
	3dcopy tmp_Intersection_GM_mask.nii.gz $mask
	rm tmp*
fi

if [ ! -f ${mask}.HEAD ]; then
	echo >&2
	echo "Could not construct $mask. Exit 3" >&2
	echo >&2; exit 3
fi


# get template
if [ ! -f vold2_mni_brain+tlrc.HEAD ]; then
	cp ${tempDir}/vold2_mni_brain+tlrc* .
fi




### --- Run MVM for each seed --- ###

for i in ${zSeeds[@]}; do

	unset subjString
	c=0; while [ $c -lt ${#subjList[@]} ]; do
		for j in ${phaseArr[@]}; do
			subjString+="${subjList[$c]} ${groupList[$c]} $j ${workDir}/${subjList[$c]}/ses-${j}/\"ZTrans_${i}+tlrc\" "			###??? Researcher Input G
		done
		let c=$[$c+1]
	done


###??? Researcher Input H
cat > MVM_${i}.sh << EOF
module load r/3.6

3dMVM -prefix MVM_${i} \
-jobs 10 \
-mask $mask \
-bsVars 'Group' \
-wsVars 'Phase' \
-num_glt 6 \
-gltLabel 1 MePhase_I-C -gltCode 1 'Group: 1*I -1*C' \
-gltLabel 2 MeGroup_W-S -gltCode 2 'Phase: -1*Sleep 1*Wake' \
-gltLabel 3 I.W-S -gltCode 3 'Group: 1*I Phase: 1*Wake -1*Sleep' \
-gltLabel 4 C.W-S -gltCode 4 'Group: 1*C Phase: 1*Wake -1*Sleep' \
-gltLabel 5 S.I-C -gltCode 5 'Group: 1*I -1*C Phase: 1*Sleep' \
-gltLabel 6 W.I-C -gltCode 6 'Group: 1*I -1*C Phase: 1*Wake' \
-dataTable \
Subj Group Phase InputFile \
$subjString
EOF

	if [ ! -f MVM_${i}+tlrc.HEAD ]; then
		source MVM_${i}.sh
	fi

	if [ ! -f MVM_${i}+tlrc.HEAD ]; then
		echo >&2
		echo "MVM faied on $i. Exit 4" >&2
		echo >&2; exit 4
	fi
done




### --- Run Monte Carlo Simulations --- ###

if [ ! -s ${outDir}/MC_output.txt ]; then

	printRaw=${outDir}/MC_acf.txt
	> $printRaw

	for i in ${subjList[@]}; do
		for j in ${phaseArr[@]}; do
			tail -n 1 ${workDir}/${i}/ses-${j}/blur_errts.1D >> $printRaw
		done
	done

	xA=`awk '{ total += $1 } END { print total/NR }' $printRaw`
	xB=`awk '{ total += $2 } END { print total/NR }' $printRaw`
	xC=`awk '{ total += $3 } END { print total/NR }' $printRaw`

	3dClustSim -mask $mask -LOTS -iter 10000 -acf $xA $xB $xC > MC_output.txt

	if [ ! -s MC_output.txt ]; then
		echo >&2
		echo "3dClustSim failed. Exit 5" >$2
		echo >&2; exit 5
	fi
fi
