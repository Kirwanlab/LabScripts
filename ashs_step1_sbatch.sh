#!/bin/bash


# Written by Nathan Muncy on 2/1/18


#SBATCH --time=05:00:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=16gb   # memory per CPU core
#SBATCH -J "ashs"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE



### Notes
#
# 1) update $tempDir for your environment
#
# 2) notice that this script calls for ashs_patch.sh - put ashs_patch.sh in $ASHS_ROOT
#
# 3) make sure the template points to the "final" dir, which contains ashs_system_config.sh


tempDir=/fslhome/<somehwere>/ashs_fast_upenn/ashs_atlas_upennpmc_20161128
dataDir=${1}/${2}

ashs_patch.sh \
-I $2 \
-a $tempDir \
-g ${dataDir}/struct_t1.nii.gz \
-f ${dataDir}/struct_t2.nii.gz \
-w ${dataDir}
