#!/bin/bash

#SBATCH --time=30:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "PPI2"   # job name

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
# 1) Will do REML deconvolution (GLS), post files, and print out info
#
# 2) Can do a variable number of deconvolutions for each phase of the experiment
#		generates the scripts, and then runs them - turn off $runDecons for simple script generation
#
# 3) Assumes timing files exist in derivatives/sub-123/timing_files
#		so a wildcard can catch only timing 1D files
#
# 4) deconNum = number of planned decons per PHASE of experimet.
#		Number of positions corresponds to step1 $phaseArr
#		Enter number of desired decons for each pahse of experiment
#			e.g. deconNum=(1 2 1) means one decon from first phase, two from second, etc.
#		Note - A value must be entered for each phase, even if no decon is desired
#			e.g. deconNum=(0 1 0) for one only decon from second phase, and no others
#
# 5) 04/17/19 - This is has been adapted from our Task scripts to be used for the
#		AutismOlfactory PPI pipeline - some things are hardcoded in the GenDecon Function!
#		-nate



subj=$1


### --- Experimenter input --- ###
#
# Change parameters for your study in this section.

parDir=~/compute/AutismOlfactory							  			# parent dir, where derivatives is located
workDir=${parDir}/derivatives/$subj

txtFile=1														# whether timing files are in txt format (1) or 1D (0)
txtTime=0														# if txt file has block duration (1:3) for pmBLOCK (1=on)
runDecons=1														# toggle for running reml scripts and post hoc (1=on) or just writing scripts (0)

deconNum=(3)													# See Note 4 above
deconPref=(FUMC FUMvC FUvC)										# array of prefix for each planned decon (length must equal sum of $deconNum)


## For 1D timing files


##### Patched here for PPI
deconLen=(3)													# trial duration for each Phase (argument for BLOCK in deconvolution). Use when $txtFile=0 or $txtTime=0
#deconTiming=(Test1_TF_4_behVect)								# array of timing files for each planned deconvolution (length must == $deconPref)


# For txt timing files
subjHold=${subj#*-}

txtFUMC=(${subjHold}_{ENI1,RI,RP,Jit1,MASK,FBO,UBO,CA}.txt)
txtFUMvC=(${subjHold}_{ENI1,RI,RP,Jit1,CA,Odor}.txt)
txtFUvC=(${subjHold}_{ENI1,RI,RP,Jit1,CA,MASK,FUBO}.txt)


# Label beh sub-bricks, per decon
namFUMC=(ENI RI RP Jit MASK FBO UBO CA)									# "Foo" of namFoo matches a $deconPref value, one string per timing file (e.g. deconPref=(SpT1); namSpT1=(Hit CR Miss FA))
namFUMvC=(ENI RI RP Jit CA Odor)
namFUvC=(ENI RI RP Jit CA MASK FUBO)





### --- Set up --- ###
#
# Determine number of phases, and number of blocks per phase
# then set these as arrays. Set up decon arrays.
# Check the deconvolution variables are set up correctly.
# Function for writing decon script.


### determine num phases/blocks
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


### checks
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

if [ $txtFile != 1 ]; then
	if [ $numDecon != ${#deconTiming[@]} ] || [ $numDecon != ${#deconPref[@]} ]; then
		echo "" >&2
		echo "Number of planned deconvolutions != to number of timing files or prefixes. Exit 2" >&2
		echo "" >&2
		exit 2
	fi
fi

if [ $phaseLen != ${#deconNum[@]} ]; then
	echo "" >&2
	echo "Length of $phaseLen != ${#deconNum[@]}. decons for each phase are required. Exit 3" >&2
	echo "" >&2
	exit 3
fi

if [ $txtFile != 1 ]; then
	tfCount=`ls timing_files/*.01.1D | wc -l`
	if [ $tfCount == 0 ]; then
		echo "" >&2
		echo "Did not detect dir \"timing_files\" or timing_files/*01.1D in $workDir. Exit 4" >&2
		echo "" >&2
		exit 4
	fi
fi

if [ $txtTime != 1 ] && [ ${#deconLen[@]} -lt 1 ]; then
	echo >&2
	echo "A block duration argument is needed if txtTime=0. Exit 5" >&2
	echo >&2
	exit 5
fi



### Function - write deconvolution script
#
# this is a lightly adjusted version of our standard decon function

GenDecon (){

	# assign vars for readability
	if [ $txtFile == 1 ]; then
		if [ $txtTime == 1 ]; then

			local h_phase=$1
			local h_block=$2
		    local h_input=$3
		    local h_out=$4
		    local h_len=$5

		    shift 5
		    local h_arr=( "$@" )
		    local nam=(${h_arr[@]:0:$h_len})
		    local txt=(${h_arr[@]:$h_len})

		else
			local h_phase=$1
			local h_block=$2
		    local h_input=$3
		    local h_out=$4
		    local h_len=$5
		    local h_trlen=$6

		    shift 6
		    local h_arr=( "$@" )
		    local nam=(${h_arr[@]:0:$h_len})
		    local txt=(${h_arr[@]:$h_len})
		fi
	else
		local h_phase=$1
		local h_block=$2
		local h_tfile=$3
		local h_trlen=$4
	    local h_input=$5
	    local h_out=$6

	    shift 6
	    local nam=( "$@" )
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


	## build behavior list
	unset stimBeh


	# if txt files supplied
	if [ $txtFile == 1 ]; then
		if [ $txtTime == 1 ]; then
			cc=0; while [ $cc -lt ${#txt[@]} ]; do
				stimBeh+="-stim_times_AM1 $x timing_files/${txt[$cc]} \"dmBLOCK(1)\" -stim_label $x beh_${nam[$cc]} "
				let x=$[$x+1]
				let cc=$[$cc+1]
			done
		else
			cc=0; while [ $cc -lt ${#txt[@]} ]; do

				###### Patched here for AO PPI study
				if [ ${txt[$cc]} == ${subjHold}_ENI1.txt ]; then
					stimBeh+="-stim_times $x timing_files/${txt[$cc]} \"BLOCK(0.5,1)\" -stim_label $x beh_${nam[$cc]} "
				elif [ ${txt[$cc]} == ${subjHold}_RI.txt ] || [ ${txt[$cc]} == ${subjHold}_RP.txt ]; then
					stimBeh+="-stim_times $x timing_files/${txt[$cc]} \"BLOCK(6,1)\" -stim_label $x beh_${nam[$cc]} "
				else
					stimBeh+="-stim_times_AM2 $x timing_files/${txt[$cc]} \"dmBLOCK(1)\" -stim_label $x beh_${nam[$cc]} "
				fi

				#stimBeh+="-stim_times $x timing_files/${txt[$cc]} \"BLOCK(${h_trlen},1)\" -stim_label $x beh_${nam[$cc]} "
				let x=$[$x+1]
				let cc=$[$cc+1]
			done
		fi

	# if 1D files supplies
	else
	    tBeh=`ls timing_files/${h_tfile}* | wc -l`
	    for ((t=1; t<=$tBeh; t++)); do
	        stimBeh+="-stim_times $x timing_files/${h_tfile}.0${t}.1D \"BLOCK(${h_trlen},1)\" -stim_label $x beh_${nam[$(($t-1))]} "
	        let x=$[$x+1]
	    done
    fi


	# num_stimts
    h_nstim=$(($x-1))


	# write script
    echo "3dDeconvolve \
    -x1D_stop \
    -input $h_input \
    -censor censor_${h_phase}_combined.1D \
    -polort A -float \
    -num_stimts $h_nstim \
    $stimBase \
    $stimBeh \
    -jobs 6 \
    -x1D X.${h_out}.xmat.1D \
    -xjpeg X.${h_out}.jpg \
    -x1D_uncensored X.${h_out}.nocensor.xmat.1D \
    -bucket ${h_out}_stats \
    -cbucket ${h_out}_cbucket \
    -errts ${h_out}_errts" > ${h_out}_deconv.sh
}




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
# A deconvolution script (foo_deconv.sh) is generated and ran for
# each planned deconvolution.

c=0; count=0; while [ $c -lt $phaseLen ]; do

	phase=${phaseArr[$c]}

	# create input list
	unset input
	for j in run-*${phase}_scale+tlrc.HEAD; do
		input+="${j%.*} "
	done

	# for each planned decon
	numD=${deconNum[$c]}
	for(( i=1; i<=$numD; i++)); do

		out=${deconPref[$count]}

		# write script
		if [ $txtFile == 1 ]; then

			holdName=($(eval echo \${nam${out}[@]}))
			holdTxt=$(eval echo \${txt${out}[@]})

			if [ $txtTime == 1 ]; then
				GenDecon $phase ${blockArr[$c]} "$input" $out ${#holdName[@]} ${holdName[@]} $holdTxt
			else
				GenDecon $phase ${blockArr[$c]} "$input" $out ${#holdName[@]} ${deconLen[$c]} ${holdName[@]} $holdTxt
			fi
		else

			holdName=($(eval echo \${nam${out}[@]}))
			GenDecon $phase ${blockArr[$c]} ${deconTiming[$count]} ${deconLen[$c]} "$input" $out ${#holdName[@]}
		fi

		# run script
		if [ -f ${out}_stats.REML_cmd ]; then
			rm ${out}_stats.REML_cmd
		fi
		source ${out}_deconv.sh

		count=$(($count+1))
	done

	let c=$[$c+1]
done




#### --- REML and Post Calcs --- ###
#
# REML deconvolution (GLS) is run, excluding WM signal.
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
	if [ ! -f ${phase}_WMe_rall+tlrc.HEAD ]; then

		3dTcat -prefix tmp_allRuns_${phase} run-*${phase}_volreg_clean+tlrc.HEAD
		3dcalc -a tmp_allRuns_${phase}+tlrc -b final_mask_WM_eroded+tlrc -expr "a*bool(b)" -datum float -prefix tmp_allRuns_${phase}_WMe
		3dmerge -1blur_fwhm 20 -doall -prefix ${phase}_WMe_rall tmp_allRuns_${phase}_WMe+tlrc
	fi


	for j in ${regArr[@]}; do
		if [ $runDecons == 1 ]; then


			# REML
			if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
				tcsh -x ${j}_stats.REML_cmd -dsort ${phase}_WMe_rall+tlrc
			fi


			# kill if REMl failed
			if [ ! -f ${j}_stats_REML+tlrc.HEAD ]; then
				echo "" >&2
				echo "REML failed on $j probably due to behavioral timing file issues ... Exit 5" >&2
				echo "" >&2
				exit 5
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
# is the number of TRs censored, which steps3-5 use. This involves
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
			dset=${k}_stats_REML+tlrc
			cp X.${k}.xmat.1D X.xmat.1D
			3dcopy ${k}_errts_REML+tlrc errts.${subj}+tlrc


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
