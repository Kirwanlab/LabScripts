#!/bin/bash


# Written by Nathan Muncy on 2/6/18


#SBATCH --time=25:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16gb   # memory per CPU core
#SBATCH -J "FSsub"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE



### Notes:
#
# 1) Ref https://surfer.nmr.mgh.harvard.edu/fswiki/HippocampalSubfields
#
# 2) Assumes proper FS directory tree is in place
#       e.g. T1 = $subjDir/mri/orig/001.mgz
#
# 3) If using a T2, point the script to it: $t2Path
#
# 4) All three options (T1, T1T2, T2) are supplied, but only one needs to be run
#
# 5) "-all" has to be run before subregion analysis can be conducted
#       if FS has already been run, this flag can be ommitted.



subj=$1
subjDir=$2
t2Path=${subjDir}/$subj

recon-all \
-all \
-s $subj \
-sd $subjDir \
-hippocampal-subfields-T1 \
-hippocampal-subfields-T1T2 ${t2Path}/T2.mgz T1T2 \
-hippocampal-subfields-T2 ${t2Path}/T2.mgz T2 \
-itkthreads 4
