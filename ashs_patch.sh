#!/bin/bash
#$ -S /bin/bash


# Two sections that have been causing the script to exit on the supercomputer have been commented out.
# Nathan Muncy 2/1/18



#######################################################################
#
#  Program:   ASHS (Automatic Segmentation of Hippocampal Subfields)
#  Module:    $Id$
#  Language:  BASH Shell Script
#  Copyright (c) 2012 Paul A. Yushkevich, University of Pennsylvania
#  
#  This file is part of ASHS
#
#  ASHS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details. 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#######################################################################

# TODO:
#   - Check that the data are in the right orientation, or handle 
#     various orientations

function usage()
{
  cat <<-USAGETEXT
		ashs_main: automatic segmentation of hippocampal subfields
		usage:
		  ashs_main [options]

		required options:
		  -a dir            Location of the atlas directory. Can be a full pathname or a
		                    relative directory name under ASHS_ROOT/data directory. 
		  -g image          Filename of 3D (g)radient echo MRI (ASHS_MPRAGE, T1w)
		  -f image          Filename of 2D focal (f)ast spin echo MRI (ASHS_TSE, T2w)
		  -w path           Working/output directory

		optional:
		  -d                Enable debugging
		  -h                Print help
		  -s integer        Run only one stage (see below); also accepts range (e.g. -s 1-3)
		  -N                No overriding of ANTS/FLIRT results. If a result from an earlier run
		                    exists, don't run ANTS/FLIRT again
		  -T                Tidy mode. Cleans up files once they are unneeded. The -N option will
		                    have no effect in tidy mode, because ANTS/FLIRT results will be erased.
		  -I string         Subject ID (for stats output). Defaults to last word of working dir.
		  -V                Display version information and exit
		  -C file           Configuration file. If not passed, uses $ASHS_ROOT/bin/ashs_config.sh
		  -Q                Use Sun Grid Engine (SGE) to schedule sub-tasks in each stage. By default,
		                    the whole ashs_main job runs in a single process. If you are doing a lot
		                    of segmentations and have SGE, it is better to run each segmentation 
		                    (ashs_main) in a separate SGE job, rather than use the -q flag. The -q flag
		                    is best for when you have only a few segmentations and want them to run fast.
		  -q OPTS           Pass in additional options to SGE's qsub. Also enables -Q option above.
      -P                Use GNU parallel to run on multiple cores on the local machine. You need to
                        have GNU parallel installed.
		  -r files          Compare segmentation results with a reference segmentation. The parameter
		                    files should consist of two nifti files in quotation marks:

		                      -r "ref_seg_left.nii.gz ref_seg_right.nii.gz"
                        
		                    The results will include overlap calculations between different
		                    stages of the segmentation and the reference segmentation. Note that the
		                    comparison takes into account the heuristic rules specified in the altas, so
		                    it is not as simple as computing dice overlaps between the reference seg
		                    and the ASHS segs.

		stages:
		  1:                fit to population template
		  2:                multi-atlas registration
		  3:                consensus segmentation using voting
		  4:                bootstrap registration
		  5:                bootstrap segmentation using voting
		  6:                segmentation Q/A
		  7:                volumes and statistics

		notes:
		  The ASHS_TSE image slice direction should be z. In other words, the dimension
		  of ASHS_TSE image should be 400x400x30 or something like that, not 400x30x400
	USAGETEXT
}

# Dereference a link - different calls on different systems
function dereflink ()
{
  if [[ $(uname) == "Darwin" ]]; then
    local SLTARG=$(readlink $1)
    if [[ $SLTARG ]]; then
      echo $SLTARG
    else
      echo $1
    fi
  else
    readlink -f $1
  fi
}

# Print usage by default
if [[ $# -lt 1 ]]; then
  echo "Try $0 -h for more information."
  exit 2
fi

# Clear the variables affected by the flags
unset ATLAS ASHS_MPRAGE ASHS_TSE ASHS_WORK STAGE_SPEC
unset ASHS_SKIP_ANTS ASHS_SKIP_RIGID ASHS_TIDY ASHS_SUBJID
unset ASHS_USE_QSUB ASHS_REFSEG_LEFT ASHS_REFSEG_RIGHT ASHS_REFSEG_LIST

# Read the options
while getopts "g:f:w:s:a:q:I:C:r:NTdhVQP" opt; do
  case $opt in

    a) ATLAS=$(dereflink $OPTARG);;
    g) ASHS_MPRAGE=$(dereflink $OPTARG);;
    f) ASHS_TSE=$(dereflink $OPTARG);;
    w) ASHS_WORK=$(dereflink $OPTARG);;
    s) STAGE_SPEC=$OPTARG;;
    N) ASHS_SKIP_ANTS=1; ASHS_SKIP_RIGID=1; ;;
    T) ASHS_TIDY=1;;
    I) ASHS_SUBJID=$OPTARG;;
    Q) ASHS_USE_QSUB=1;;
    P) ASHS_USE_PARALLEL=1;;
    q) ASHS_USE_QSUB=1; QOPTS=$OPTARG;;
    C) ASHS_CONFIG=$(dereflink $OPTARG);;
    r) ASHS_REFSEG_LIST=($(echo $OPTARG));;
    d) set -x -e;;
    h) usage; exit 0;;
    V) vers; exit 0;;
    \?) echo "Unknown option $OPTARG"; exit 2;;
    :) echo "Option $OPTARG requires an argument"; exit 2;;

  esac
done

# Check the root dir
if [[ ! $ASHS_ROOT ]]; then
  echo "Please set ASHS_ROOT to the ASHS root directory before running $0"
  exit -2
#elif [[ $ASHS_ROOT != $(dereflink $ASHS_ROOT) ]]; then
#  echo "ASHS_ROOT must point to an absolute path, not a relative path"
#  exit -2
fi

# Set the config file
if [[ ! $ASHS_CONFIG ]]; then
  ASHS_CONFIG=$ATLAS/ashs_user_config.sh
fi

# Check that parallel and qsub are not both on
if [[ $ASHS_USE_PARALLEL && $ASHS_USE_QSUB ]]; then
  echo "Cannot use SGE (-Q) and GNU Parallel (-P) at the same time"
  exit -2
fi

# Load the library. This also processes the config file
source $ASHS_ROOT/bin/ashs_lib.sh

# Check if the required parameters were passed in
echo "Atlas    : ${ATLAS?    "Directory for atlas was not specified. See $0 -h"}"
echo "T1 Image : ${ASHS_MPRAGE?   "T1-weighted MRI was not specified. See $0 -h"}"
echo "T2 Image : ${ASHS_TSE?      "T2-weighted MRI was not specified. See $0 -h"}"
echo "WorkDir  : ${ASHS_WORK?     "Working directory was not specified. See $0 -h"}"

# Handle the -r parameter
if [[ $ASHS_REFSEG_LIST ]]; then
  if [[ ${#ASHS_REFSEG_LIST[*]} -eq 2 ]]; then
    ASHS_REFSEG_LEFT=$(dereflink ${ASHS_REFSEG_LIST[0]})
    ASHS_REFSEG_RIGHT=$(dereflink ${ASHS_REFSEG_LIST[1]})
    if [[ ! -f $ASHS_REFSEG_LEFT || ! -f $ASHS_REFSEG_RIGHT ]]; then
      echo "Reference segmentation $ASHS_REFSEG_LEFT $ASHS_REFSEG_RIGHT not found"
      exit -2
    fi
  else
    echo "Wrong number of parameters to -r option"
  fi
fi

# Whether we are using QSUB        ### I commented out the exit here, -Nate
if [[ $ASHS_USE_QSUB ]]; then
  if [[ ! $SGE_ROOT ]]; then
    echo "-Q flag used, but SGE is not present."
#exit -1;
  fi
  echo "Using SGE with root $SGE_ROOT and options $QOPTS"
elif [[ $ASHS_USE_PARALLEL ]]; then
  echo "Using GNU parallel"
else
  echo "Not using SGE"
fi

# Convert the work directory to absolute path
mkdir -p ${ASHS_WORK?}
ASHS_WORK=$(cd $ASHS_WORK; pwd)
if [[ ! -d $ASHS_WORK ]]; then 
  echo "Work directory $ASHS_WORK does not exist";
fi

# Check the atlas location
if [[ -f $ATLAS/ashs_atlas_vars.sh ]]; then
  ASHS_ATLAS=$ATLAS;
elif [[ -f $ASHS_ROOT/data/$ATLAS/ashs_atlas_vars.sh ]]; then
  ASHS_ATLAS=$ASHS_ROOT/data/$ATLAS
else
  echo "Atlas directory must be specified"
  exit 2;
fi

# Check the heuristics in the atlas
if [[ -f $ASHS_ATLAS/ashs_heuristics.txt ]]; then
	ASHS_HEURISTICS=$ASHS_ATLAS/ashs_heuristics.txt
fi

# Make sure all files exist
if [[ ! $ASHS_MPRAGE || ! -f $ASHS_MPRAGE ]]; then
	echo "T1-weighted 3D gradient echo MRI (-g) must be specified"
	exit 2;
elif [[ ! $ASHS_TSE || ! -f $ASHS_TSE ]]; then
	echo "T2-weighted 2D fast spin echo MRI (-f) must be specified"
	exit 2;
elif [[ ! $ASHS_WORK ]]; then
	echo "Working/output directory must be specified"
	exit 2;
fi

# Check that the dimensions of the T2 image are right
DIMS=$(c3d $ASHS_TSE -info | cut -d ';' -f 1 | sed -e "s/.*\[//" -e "s/\].*//" -e "s/,//g")
if [[ ${DIMS[2]} > ${DIMS[0]} || ${DIMS[2]} > ${DIMS[1]} ]]; then
  echo "The T2-weighted image has wrong dimensions (fails dim[2] < min(dim[0], dim[1])"
  exit -1
fi

# Subject ID set to work dir last work
if [[ ! $ASHS_SUBJID ]]; then
  ASHS_SUBJID=$(basename $ASHS_WORK)
fi

# Create the working directory and the dump directory
mkdir -p $ASHS_WORK $ASHS_WORK/dump $ASHS_WORK/final

# Set the start and end stages
if [[ $STAGE_SPEC ]]; then
  STAGE_START=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $1}')
  STAGE_END=$(echo $STAGE_SPEC | awk -F '-' '$0 ~ /^[0-9]+-*[0-9]*$/ {print $NF}')
else
  STAGE_START=1
  STAGE_END=15
fi

if [[ ! $STAGE_END || ! $STAGE_START ]]; then
  echo "Wrong stage specification -s $STAGE_SPEC"
  exit -1;
fi

# Get the number of atlases, other information
source $ASHS_ATLAS/ashs_atlas_vars.sh

# List of sides for the array qsub commands below
SIDES="$ASHS_SIDES"

# Run the stages of the script
export ASHS_ROOT ASHS_WORK ASHS_SKIP_ANTS ASHS_SKIP_RIGID ASHS_SUBJID ASHS_CONFIG ASHS_ATLAS
export ASHS_HEURISTICS ASHS_TIDY ASHS_MPRAGE ASHS_TSE ASHS_REFSEG_LEFT ASHS_REFSEG_RIGHT QOPTS
export SIDES

# List of training atlases 
TRIDS=$(for((i = 0; i < $ASHS_ATLAS_N; i++)); do echo $(printf "%03i" $i); done)

for ((STAGE=$STAGE_START; STAGE<=$STAGE_END; STAGE++)); do

  case $STAGE in 

    1) 
    # Template matching
    echo "Running stage 1: normalize to T1 population template"
    qsubmit_sync "ashs_stg1" $ASHS_ROOT/bin/ashs_template_qsub.sh ;;

    2) 
    # Multi-atlas matching 
    echo "Running stage 2: normalize to multiple T1/T2 atlases"
    qsubmit_double_array "ashs_stg2" "$SIDES" "$TRIDS" $ASHS_ROOT/bin/ashs_multiatlas_qsub.sh ;;

    3) 
    # Voting
    echo "Running stage 3: Label Fusion"
    qsubmit_single_array "ashs_stg3" "$SIDES" $ASHS_ROOT/bin/ashs_voting_qsub.sh 0 ;;

    4)
    # Bootstrapping
    echo "Running stage 4: Bootstrap segmentation"
    qsubmit_double_array "ashs_stg4" "$SIDES" "$TRIDS" $ASHS_ROOT/bin/ashs_bootstrap_qsub.sh ;;

    5)
    # Bootstrap voting
    echo "Running stage 5: Bootstrap label fusion" 
    qsubmit_single_array "ashs_stg5" "$SIDES" $ASHS_ROOT/bin/ashs_voting_qsub.sh 1 ;;

    6)
    # Final QA
    echo "Running stage 6: Final QA"
    qsubmit_sync "ashs_stg6" $ASHS_ROOT/bin/ashs_finalqa_qsub.sh ;;
  
    7) 
    # Statistics & Volumes
    echo "Running stage 7: Statistics and Volumes"
    qsubmit_sync "ashs_stg7" $ASHS_ROOT/bin/ashs_extractstats_qsub.sh ;;

  esac  

done
