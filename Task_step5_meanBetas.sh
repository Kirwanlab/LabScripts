#!/bin/bash

#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=6gb   # memory per CPU core
#SBATCH -J "TS5"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE






# Written by Nathan Muncy on 11/28/18


###--- Notes, in no particular order
#
# 1) the script will split the cluster file into multiple masks, and pull betas from each participant.
#
# 2) assumes clusters from step4 output have been saved in Clust_$fileArr format
#		will use the comparisonString portion to keep different group analyses straight
#		comparisonString should match the decon prefix (step4 $etacList)
#
# 3) assumes that decon files exist locally (bring them back from the supercomputer)
#
# 4) Written for the output of ETAC - will pull betas for each appropriate blur from each cluster



# Variables
parDir=~/compute/STT_reml
workDir=${parDir}/derivatives										###??? Update this section
grpDir=${parDir}/Analyses/grpAnalysis
clustDir=${grpDir}/etac_clusters
outDir=${grpDir}/etac_betas
refDir=${workDir}/sub-1295											# reference file for dimensions etc


fileArr=(SpT1 SpT1pT2 T1 T1pT2 T2 T2fT1)							# decon files from which betas will be extracted - should match step4.
arrA=(1 7 1 7 1 1)													# sub-bricks corresponding to $fileArr
arrB=(4 10 4 10 4 7)
arrLen=${#arrA[@]}

blurX=({2..4})														# blur mulitpliers from previous steps




# function - search array for string
MatchString (){
	local e match="$1"
	shift
	for e; do [[ "$e" == "$match" ]] && return 0; done
	return 1
}


# check
test1=${#fileArr[@]}; test2=${#arrB[@]}
if [ $test1 != $arrLen ] || [ $test2 != $arrLen ]; then
	echo "Script set up incorreclty. Exit 1" >&2
	exit 1
fi


# determine blur
refFile=${refDir}/${fileArr[0]}_stats_REML+tlrc
gridSize=`3dinfo -dk $refFile`
int=`printf "%.0f" $gridSize`

c=0; for i in ${blurX[@]}; do
	hold="$(($int * $i))"
	blurArr[$c]=$hold
	let c=$[$c+1]
done
blurLen=${#blurArr[@]}



### make clusters, tables
mkdir $outDir $clustDir
cd ${grpDir}/etac_indiv

for i in FINAL_*allP*.HEAD; do

	tmp1=${i%_ET*}; pref=${tmp1##*_}
	tmp2=${i#*_}; blur=${tmp2%%_*}

	if [ ! -f ${clustDir}/Clust_${pref}_${blur}_table.txt ]; then

		3dclust -1Dformat -nosum -1dindex 0 \
		-1tindex 0 -2thresh -0.5 0.5 -dxyz=1 \
		-savemask Clust_${pref}_${blur}_mask \
		1.01 5 $i > Clust_${pref}_${blur}_table.txt
	fi
done

# organize
mv Clust* $clustDir



### pull mean betas for e/cluster from e/comparison from e/subject
cd $clustDir

c=0; while [ $c -lt $arrLen ]; do

	hold=${fileArr[$c]}
	betas=${arrA[$c]},${arrB[$c]}

	# make subj list
	unset subjHold
	arrRem=(`cat ${grpDir}/info_rmSubj_${hold}.txt`)
	for i in ${workDir}/s*; do
		subj=${i##*\/}
		MatchString "$subj" "${arrRem[@]}"
		if [ $? == 1 ]; then
			subjHold+="$subj "
		fi
	done
	subjList=(${subjHold})


	# split clust masks
	d=0; while [ $d -lt $blurLen ]; do

		blurInt=${blurArr[$d]}

		if [ -f Clust_${hold}_b${blurInt}_mask+tlrc.HEAD ]; then
			if [ ! -f Clust_${hold}_b${blurInt}_c1+tlrc.HEAD ]; then

				3dcopy Clust_${hold}_b${blurInt}_mask+tlrc ${hold}_b${blurInt}.nii.gz
				num=`3dinfo Clust_${hold}_b${blurInt}_mask+tlrc | grep "At sub-brick #0 '#0' datum type is short" | sed 's/[^0-9]*//g' | sed 's/^...//'`

				for (( j=1; j<=$num; j++ )); do
					if [ ! -f Clust_${hold}_b${blurInt}_c${j}+tlrc.HEAD ]; then

						c3d ${hold}_b${blurInt}.nii.gz -thresh $j $j 1 0 -o ${hold}_b${blurInt}_${j}.nii.gz
						3dcopy ${hold}_b${blurInt}_${j}.nii.gz Clust_${hold}_b${blurInt}_c${j}+tlrc
					fi
				done
				rm *.nii.gz
			fi


			# pull betas
			for i in Clust_${hold}_b${blurInt}_c*+tlrc.HEAD; do

				tmp=${i##*_}
				cnum=${tmp%+*}
				print=${outDir}/Betas_${hold}_b${blurInt}_${cnum}.txt
				> $print

				for j in ${subjList[@]}; do

					subjDir=${workDir}/${j}
					decon=stats_REML

					# blur
					if [ ! -f ${subjDir}/${hold}_${decon}_blur${blurInt}+tlrc.HEAD ]; then
						3dmerge -prefix ${subjDir}/${hold}_${decon}_blur${blurInt} -1blur_fwhm $blurInt -doall ${subjDir}/${hold}_${decon}+tlrc
					fi

					file=${subjDir}/${hold}_${decon}_blur${blurInt}+tlrc
					stats=`3dROIstats -mask $i "${file}[${betas}]"`
					echo "$j $stats" >> $print
				done
			done
		fi
		let d=$[$d+1]
	done
	let c=$[$c+1]
done
