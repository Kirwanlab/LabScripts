#!/bin/bash

#SBATCH --time=5:00:00   # walltime
#SBATCH --ntasks=2   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4gb   # memory per CPU core
#SBATCH -J "RS3"   # job name

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
# This script geneartes seeds, pulls correlation matrix,
#	and Fischer Z-transforms the matrices
#
# Major steps following the comment ### --- Foo --- ###,
# minor steps are also annotated
#
# Sections ###??? require research attention
#
# Researcher Input A) Array of desired JLF priorss. Length must equal length of roiName
# Researcher Input B) Array of desired ROI names
# Researcher Input C) Array of coordinates for seed construction. Length must equal length of coordName
# Researcher Input D) Array of desired seed names for coordinates
# Researcher Input E) Desired seed size (radius, integer)



subj=$1
phase=$2




### --- Set Up --- ###

workDir=~/compute/MNI_NREM/MRI_processed/derivatives/${subj}/ses-$phase 	# Location of subject-session directory
tempDir=~/bin/Templates/vold2_mni    				# Location of template/priors
jlfDir=${tempDir}/priors_JLF						# Priors dir


## Flags
jlfFlag=0											# if JLF priors/ROIs are being used as seeds (1=on)
corFlag=1											# if coordinates are being used for seeds (1=on)


# JLF seeds
roiList=(0053 0054)									###??? Research Input A
roiName=(RHC RAmg)									###??? Research Input B
roiLen=${#roiList[@]}


# Coordinate seeds
coordList=("5 -55 25" "-40 22 38")					###??? Research Input C
coordName=(rPCC lMFG)								###??? Research Input D
seedSize=5      									###??? Research Input E
coordLen=${#coordList[@]}




### --- Checks --- ###

if [ $jlfFlag == 1 ]; then
	if [ ${#roiList[@]} == 0 ]; then
		echo >&2
		echo "At least one JLF seed needed with jlfFlag=1. Exit 1" >&2
		echo >&2; exit 1
	fi
	if [ ${#roiList[@]} != ${#roiName[@]} ]; then
		echo >&2
		echo "Length of roiList and roiName are not equal. Exit 2" >&2
		echo >&2; exit 2
	fi
fi

if [ $corFlag == 1 ]; then
	if [ ${#coordList[@]} == 0 ]; then
		echo >&2
		echo "At least one coordinate seed needed with corFlag=1. Exit 3" >&2
		echo >&2; exit 3
	fi
	if [ ${#coordList[@]} != ${#coordName[@]} ]; then
		echo >&2
		echo "Length of coordList and coordName are not equal. Exit 4" >&2
		echo >&2; exit 4
	fi
fi




### --- Seed Timeseries --- ###
#
# Copy or construct time seed, make sure
# seed is in functional resoultion, and make
# correlation file

cd $workDir
RSfile=errts_fanaticor+tlrc

## For JLF
if [ $jlfFlag == 1 ]; then
    c=0; while [ $c -lt $roiLen ]; do

        seedName=Seed_${roiName[$c]}

        if [ ! -s FINAL_${seedName}+tlrc.HEAD ]; then

        	# get prior, resample, binarize
            jlfPrior=${jlfDir}/label_${roiList[$c]}.nii.gz

            3dcopy $jlfPrior tmp_${seedName}+tlrc
            3dfractionize -template $RSfile -input tmp_${seedName}+tlrc -prefix tmp_rs_${seedName}
            3dcalc -a tmp_rs_${seedName}+tlrc -prefix $seedName -expr "step(a)" && rm tmp*

            # extract correlation
            3dROIstats -quiet -mask ${seedName}+tlrc $RSfile > ${seedName}_TimeSeries.1D
            3dTcorr1D -mask full_mask+tlrc -prefix FINAL_${seedName} $RSfile ${seedName}_TimeSeries.1D
        fi
    let c=$[$c+1]
    done
fi


## For coordinate
if [ $corFlag == 1 ]; then
    c=0; while [ $c -lt $coordLen ]; do

        seedName=Seed_${coordName[$c]}

        if [ ! -s FINAL_${seedName}+tlrc.HEAD ]; then

        	# create seed, extract correlation
            echo ${coordList[$c]} > ${seedName}.txt
            3dUndump -prefix $seedName -master $RSfile -srad $seedSize -xyz ${seedName}.txt
            3dROIstats -quiet -mask ${seedName}+tlrc $RSfile > ${seedName}_TimeSeries.1D
            3dTcorr1D -mask full_mask+tlrc -prefix FINAL_${seedName} $RSfile ${seedName}_TimeSeries.1D
        fi
    let c=$[$c+1]
    done
fi




### --- Z-Transform --- ###

for i in FINAL*HEAD; do

	scan=${i%.*}
	tmp=${i%+*}
	seed=${tmp##*_}

	if [ ! -f ZTrans_${seed}+tlrc.HEAD ]; then
		3dcalc -a $scan -expr 'log((1+a)/(1-a))/2' -prefix ZTrans_$seed
	fi
done