#!/bin/bash

#SBATCH --time=30:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "sttR2"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE




# Written by Nathan Muncy on 10/24/18


### --- Notes
#
# 1) Will do deconvolutions (OLS, GLS), post files, and print out info
#
# 2) Can do a variable number of deconvolutions for each phase of the experiment
#		generates the scripts, and then runs them - turn off $runDecons for simple script generation
#
# 3) Assumes timing files exist in derivatives/sub-123/timing_files, and are in .1D format
#
# 4) Only the variables at the beginning of the script need to be updated (this is NEW and COOL)




subj=$1
testMode=$2


														###??? update these variables/arrays - NEW: only this section needs to be updated!
parDir=~/compute/STT_reml					  			# parent dir, where derivatives is located
workDir=${parDir}/derivatives/$subj

doREML=1												# conduct GLS decon (1=on), runDecons must be 1
runDecons=1												# toggle for running decon/reml scripts and post hoc (1=on)

deconNum=(2 2 2)										# number of planned decons per PHASE, corresponds to $phaseArr from step1 (STUDY TEST1 TEST2). E.g. (2 1 1) = the first phase will have 2 deconvolutions, and the second and third phase will both have 1 respectively
deconLen=(3 3 4.5)										# trial duration for each Phase (argument for BLOCK in deconvolution)
deconPref=(SpT1 SpT1pT2 T1 T1pT2 T2 T2fT1)				# array of prefix for each planned decon (length must equal sum of $deconNum)

deconTiming=(
Study_pred_Test1_TF_4_behVect
Study_pred_Test1_pred_Test2_TF_4_behVect
Test1_TF_4_behVect
Test1_pred_Test2_TF_4_behVect
Test2_TF_4_behVect
Test2_ref_Test1_TF_4_behVect
)														# array of timing files for each planned deconvolution (length must == $deconPref)



### --- Not Lazy Section - for those who really like their data well organized
#
# This section is for setting the behavioral sub-brick labels (beh_foo, beh_bar). If not
# used (NotLazy=0), then the labels are beh_1, beh_2, etc.
#
# In order to work, the arrays below must have a title that matches $deconPref (e.g. arrSpT1 <- deconPref=(SpT1))
#
# Also, the length of the array must be equal to the number of .1D files for the deconvolution
#
# The number of arrays should equal the length of $deconTiming

NotLazy=1												# 1=on

arrSpT1=(R-Hit R-FA R-CR R-Miss NR)						# arrFoo matches a $deconPref value, one string per .1D file (e.g. arrSpT1=(Hit CR Miss FA))
arrSpT1pT2=(R-HH R-HF R-FH R-FF R-MCH R-MCF NR)
arrT1=(T1Hit T1FA T1CR T1Miss NR)
arrT1pT2=(HH HF FH FF MH MF CH CF NR)
arrT2=(T2Hit T2FA NR)
arrT2fT1=(H1H2 H1F2 F1H2 F1F2 M1H2 M1F2 C1H2 C1F2 NR)




### --- Set up --- ###
#
# Determine number of phases, and number of blocks per phase
# then set these as arrays. Set up decon arrays.
# Check the deconvolution variables are set up correctly


# determine num phases/blocks
cd $workDir

> tmp.txt
for i in run*scale+tlrc.HEAD; do

	tmp=${i%_*}
	run=${i%%_*}
	phase=${tmp#*_}

	echo -e "$run \t $phase" >> tmp.txt
done

awk -F '\t' '{print $2}' tmp.txt | sort | uniq -c > phase_list.txt
rm tmp.txt


blockArr=(`cat phase_list.txt | awk '{print $1}'`)
phaseArr=(`cat phase_list.txt | awk '{print $2}'`)
phaseLen=${#phaseArr[@]}


# checks
if [ ! ${#phaseArr[@]} -gt 0 ]; then
	echo "" >&2
	echo "Problem determinig number of phases and blocks per phase. Check step1 setup. Exit 1" >&2
	echo "" >&2
	exit 1
fi

unset numDecon
for i in ${deconNum[@]}; do
	numDecon=$(( $numDecon + $i ))
done

if [ $numDecon != ${#deconTiming[@]} ] || [ $numDecon != ${#deconPref[@]} ]; then
	echo "" >&2
	echo "Number of planned deconvolutions != to number of timing files or prefixes. Exit 2" >&2
	echo "" >&2
	exit 2
fi

if [ $phaseLen != ${#deconNum[@]} ]; then
	echo "" >&2
	echo "Length of $phaseLen != ${#deconNum[@]}. decons for each phase are required. Exit 3" >&2
	echo "" >&2
	exit 3
fi

tfCount=`ls timing_files/*.01.1D | wc -l`
if [ $tfCount == 0 ]; then
	echo "" >&2
	echo "Did not detect dir \"timing_files\" or timing_files/*01.1D in $workDir. Exit 4" >&2
	echo "" >&2
	exit 4
fi




### --- Motion --- ###
#
# motion and censor files are constructed. Multiple motion files
# include mean and derivative of motion.


c=0; while [ $c -lt $phaseLen ]; do

	phase=${phaseArr[$c]}
	nruns=${blockArr[$c]}

	cat dfile.run-*${phase}.1D > dfile_rall_${phase}.1D

	if [ ! -s censor_${phase}_combined.1D ]; then

		# files: de-meaned, motion params (per phase)
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -demean -write motion_demean_${phase}.1D
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -derivative -demean -write motion_deriv_${phase}.1D
		1d_tool.py -infile motion_demean_${phase}.1D -set_nruns $nruns -split_into_pad_runs mot_demean_${phase}
		1d_tool.py -infile dfile_rall_${phase}.1D -set_nruns $nruns -show_censor_count -censor_prev_TR -censor_motion 0.3 motion_${phase}


		# determine censor
		cat out.cen.run-*${phase}.1D > outcount_censor_${phase}.1D
		1deval -a motion_${phase}_censor.1D -b outcount_censor_${phase}.1D -expr "a*b" > censor_${phase}_combined.1D
	fi

	let c=$[$c+1]
done




### --- Deconvolve --- ###
#
# A deconvolution (OLS) script is written (for review) and then run.
# Currently, behaviors are called beh_1, beh_2 etc due to practical
# constraints. I'll probably fix this in the future. A deconv script
# is written for each run of the experiment, and for different
# phases a conditional allows control over things like timing
# files and different block times. Timing files are not constructed.
# A REML (GLS) script is written, for future use.


# Function - write deconvolution script
GenDecon (){

	# assign vars for readability
	local h_phase=$1
	local h_block=$2
	local h_tfile=$3
	local h_trlen=$4
    local h_input=$5
    local h_out=$6

    if [ $NotLazy == 1 ]; then
	    for aa in {1..6}; do
		    shift
	    done
	    local h_arr=( "$@" )
	fi


	# build motion list
	unset stimBase
    x=1

    for ((r=1; r<=${h_block}; r++)); do
        for ((b=0; b<=5; b++)); do

            stimBase+="-stim_file $x mot_demean_${h_phase}.r0${r}.1D'[$b]' -stim_base $x -stim_label $x mot_$x "
            let x=$[$x+1]
        done
    done


	# build behavior list
	unset stimBeh
    tBeh=`ls timing_files/${h_tfile}* | wc -l`

	if [ $NotLazy == 1 ]; then
	    for ((t=1; t<=$tBeh; t++)); do
	        stimBeh+="-stim_times $x timing_files/${h_tfile}.0${t}.1D \"BLOCK(${h_trlen},1)\" -stim_label $x beh_${h_arr[$(($t-1))]} "
	        let x=$[$x+1]
	    done
	else
	    for ((t=1; t<=$tBeh; t++)); do
	        stimBeh+="-stim_times $x timing_files/${h_tfile}.0${t}.1D \"BLOCK(${h_trlen},1)\" -stim_label $x beh_$t "
	        let x=$[$x+1]
	    done
    fi


	# num_stimts
    h_nstim=$(($x-1))


	# write script
    echo "3dDeconvolve \
    -input $h_input \
    -censor censor_${h_phase}_combined.1D \
    -polort A -float \
    -num_stimts $h_nstim \
    $stimBase \
    $stimBeh \
    -jobs 6 \
    -bout -fout -tout \
    -x1D X.${h_out}.xmat.1D \
    -xjpeg X.${h_out}.jpg \
    -x1D_uncensored X.${h_out}.nocensor.xmat.1D \
    -errts ${h_out}_errts \
    -bucket ${h_out}_stats" > ${h_out}_deconv.sh
}


### write, run deconvolution scripts
c=0; count=0; while [ $c -lt $phaseLen ]; do

	phase=${phaseArr[$c]}
	trialLen=${deconLen[$c]}

	# create input list
	unset input
	for j in run-*${phase}_scale+tlrc.HEAD; do
		input+="${j%.*} "
	done


	# for each planned decons
	numD=${deconNum[$c]}
	for(( i=1; i<=$numD; i++)); do

		# pull timing, prefix
		out=${deconPref[$count]}
		tfile=${deconTiming[$count]}


		# write script
		if [ $NotLazy == 1 ]; then
			hold=$(eval echo \${arr${out}[@]})
			GenDecon $phase ${blockArr[$c]} $tfile $trialLen "$input" $out $hold
		else
			GenDecon $phase ${blockArr[$c]} $tfile $trialLen "$input" $out
		fi
		count=$(($count+1))

		# run script
		if [ $runDecons == 1 ]; then
			if [ ! -f ${out}_stats+tlrc.HEAD ]; then
				source ${out}_deconv.sh
			fi
		fi
	done
	let c=$[$c+1]
done




#### --- REML and Post Calcs --- ###
#
# REML deconvolution (GLS) is run, excluding WM signal. REML will
# probably become standard soon, so I'll get this working at some point.
# Global SNR and corr are calculated.


c=0; count=0; while [ $c -lt $phaseLen ]; do


	# loop through number of planned decons, set arr
	phase=${phaseArr[$c]}

	numD=${deconNum[$c]}
	x=0; for((i=1; i<=$numD; i++)); do

		regArr[$x]=${deconPref[$count]}

		let x=$[$x+1]
		let count=$[$count+1]
	done


	# all runs signal
	countS=`1d_tool.py -infile censor_${phase}_combined.1D -show_trs_uncensored encoded`

	if [ ! -f ${regArr[0]}_TSNR+tlrc.HEAD ]; then
		3dTcat -prefix tmp_${phase}_all_runs run-*${phase}_scale+tlrc.HEAD
		3dTstat -mean -prefix tmp_${phase}_allSignal tmp_${phase}_all_runs+tlrc"[${countS}]"
	fi


	# timeseries of eroded WM
	if [ $doREML == 1 ]; then
		if [ ! -f ${phase}_WMe_rall+tlrc.HEAD ]; then

			3dTcat -prefix tmp_allRuns_${phase} run-*${phase}_volreg_clean+tlrc.HEAD
			3dcalc -a tmp_allRuns_${phase}+tlrc -b final_mask_WM_eroded+tlrc -expr "a*bool(b)" -datum float -prefix tmp_allRuns_${phase}_WMe
			3dmerge -1blur_fwhm 20 -doall -prefix ${phase}_WMe_rall tmp_allRuns_${phase}_WMe+tlrc
		fi
	fi


	for j in ${regArr[@]}; do

		# kill if decon failed
		if [ $runDecons == 1 ]; then
			if [ ! -f ${j}_stats+tlrc.HEAD ]; then
				echo "" >&2
				echo "Decon failed on $j ... Exit 5" >&2
				echo "" >&2
				exit 5
			fi
		fi


		if [ $runDecons == 1 ]; then
			if [ $doREML == 1 ]; then

				# REML
				if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
					tcsh -x ${j}_stats.REML_cmd -dsort ${phase}_WMe_rall+tlrc
				fi


				# kill if REMl failed
				if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
					echo "" >&2
					echo "REML failed on $j ... Exit 6" >&2
					echo "" >&2
					exit 6
				fi


				# calc SNR, corr
				if [ ! -f ${j}_TSNR+tlrc.HEAD ]; then

					3dTstat -stdev -prefix tmp_${j}_allNoise ${j}_errts_REML+tlrc"[${countS}]"

					3dcalc -a tmp_${phase}_allSignal+tlrc \
					-b tmp_${j}_allNoise+tlrc \
					-c full_mask+tlrc \
					-expr 'c*a/b' -prefix ${j}_TSNR

					3dTnorm -norm2 -prefix tmp_${j}_errts_unit ${j}_errts_REML+tlrc
					3dmaskave -quiet -mask full_mask+tlrc tmp_${j}_errts_unit+tlrc > ${j}_gmean_errts_unit.1D
					3dcalc -a tmp_${j}_errts_unit+tlrc -b ${j}_gmean_errts_unit.1D -expr 'a*b' -prefix tmp_${j}_DP
					3dTstat -sum -prefix ${j}_corr_brain tmp_${j}_DP+tlrc
				fi
			fi
		fi


		# detect pairwise cor
		1d_tool.py -show_cormat_warnings -infile X.${j}.xmat.1D | tee out.${j}.cormat_warn.txt
	done
	let c=$[$c+1]
done



for i in ${deconPref[@]}; do

	# sum of regressors, stim only x-matrix
	if [ ! -s X.${i}.stim.xmat.1D ]; then

		reg_cols=`1d_tool.py -infile X.${i}.nocensor.xmat.1D -show_indices_interest`
		3dTstat -sum -prefix ${i}_sum_ideal.1D X.${i}.nocensor.xmat.1D"[$reg_cols]"
		1dcat X.${i}.nocensor.xmat.1D"[$reg_cols]" > X.${i}.stim.xmat.1D
	fi
done




#### --- Print out info, Clean --- ###
#
# Print out information about the data - of particulatr interest
# is the number of TRs censored, which steps3/4 use. This involves
# producing the needed files, generating a set of review scripts,
# and then running my favorite. Many intermediates get removed.


# organize files for what gen*py needs
3dcopy full_mask+tlrc full_mask.${subj}+tlrc

if [ $runDecons == 1 ]; then
	for i in ${phaseArr[@]}; do

		cat outcount.run-*${i}.1D > outcount_all_${i}.1D

		c=1; for j in run-*${i}*+orig.HEAD; do

			prefix=${j%+*}
			3dcopy ${j%.*} pb00.${subj}.r0${c}.tcat
			3dcopy ${prefix}_volreg_clean+tlrc pb02.${subj}.r0${c}.volreg

			let c=$[$c+1]
		done


		for k in ${deconPref[@]}; do

			# a touch more organization (gen*py is very needy)
			dset=${k}_stats+tlrc
			cp X.${k}.xmat.1D X.xmat.1D
			3dcopy ${k}_errts+tlrc errts.${subj}+tlrc


			# generate script
			gen_ss_review_scripts.py \
			-subj ${subj} \
			-rm_trs 0 \
			-motion_dset dfile_rall_${i}.1D \
			-outlier_dset outcount_all_${i}.1D \
			-enorm_dset  motion_${i}_enorm.1D \
			-mot_limit 0.3 \
			-out_limit 0.1 \
			-xmat_regress X.${k}.xmat.1D \
			-xmat_uncensored X.${k}.nocensor.xmat.1D \
			-stats_dset ${dset} \
			-final_anat final_anat+tlrc \
			-final_view tlrc \
			-exit0


			# run script - write an output for e/analysis
			./\@ss_review_basic | tee out_summary_${k}.txt


			# clean
			rm errts.*
			rm X.xmat.1D
			rm pb0*
			rm *ss_review*
		done
	done
fi



# clean
if [ $testMode == 1 ]; then
	rm tmp_*
	rm -r a*
	rm final_mask_{CSF,GM}*
	rm *corr_brain*
	rm *gmean_errts*
	rm *volreg*
	rm Temp*
	rm *WMe_rall*
	rm full_mask.*
fi
