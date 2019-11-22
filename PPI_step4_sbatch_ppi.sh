#!/bin/bash

#SBATCH --time=20:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "PPI4"   # job name

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE






### Notes:
#
# This script will extract the behavioral timeseries, and generate
# correlation matrices for each seed.
#
# Somethings are specific to the AutismOlfactory study, like everything
# between lines 287-418, so update accordingly.




### --- Set up --- ###
#
# This section will orient the rest of the script. Aside from the
# contrast matrices section (parts of step1,2) this is the only
# section that needs to be updated. Future versions may be written
# more robustly to minimize the need for research input.


# General Variables
workDir=~/compute/AutismOlfactory
tempDir=${workDir}/Template
scriptDir=${workDir}/Scripts
timingDir=${workDir}/derivatives/TimingFiles


# Subject Variables
subj=$1
string=${subj#*-}
ppiDir=${workDir}/derivatives/${subj}


## Deconvolution variables												# Prefix of decons run in step2
deconList=(FUMC FUMvC FUvC)

# For txt timing files													# Same as in step2
txtFUMC=(${string}_{ENI1,RI,RP,Jit1,MASK,FBO,UBO,CA}.txt)
txtFUMvC=(${string}_{ENI1,RI,RP,Jit1,CA,Odor}.txt)
txtFUvC=(${string}_{ENI1,RI,RP,Jit1,CA,MASK,FUBO}.txt)

# Label beh sub-bricks, per decon										# Same as in step2
namFUMC=(ENI RI RP Jit MASK FBO UBO CA)
namFUMvC=(ENI RI RP Jit CA Odor)
namFUvC=(ENI RI RP Jit CA MASK FUBO)


## PPI variables														# All desirable sub-brick labels from step2 output
behInterest=(MASK FBO UBO CA Odor)

arrFUMC=(Mask FBO UBO CA)												# Behaviors to be extracted from each TS, per decon
arrFUMvC=(CA Odor)
arrFUvC=(CA Mask FUBO)

seedCoord=("-29 1 -18" "27 2 -21" "20 -3 -14" "-23 -4 -13" "46 44 23" "20 -81 41" "-32 -76 39" "-42 -75 36" "-66 -12 25" "57 -54 -12")												# Seed coordinates
seedName=(LPF RPF 1 2 5 7 8a 8b 9 10)									# Seed prefix
seedLen=${#seedCoord[@]}




### --- Functions --- ###

# Search array for string
MatchString () {

	local e match="$1"

	shift
	for e; do
		[[ "$e" == "$match" ]] && return 0
	done
	return 1
}


# Write deconvolution script
#
# This is a stripped version of the function in PPI_step2.
# It has been adjusted for use with the PPI, and written
# Specifically for the AutismOlfactory study

GenDecon (){

	# assign vars for readability
	local h_phase=$1
	local h_block=$2
    local h_input=$3
    local h_out=$4
    local h_len=$5

    shift 5
    local h_arr=( "$@" )
    local nam=(${h_arr[@]:0:$h_len})
    local txt=(${h_arr[@]:$h_len})


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

	cc=0; while [ $cc -lt ${#txt[@]} ]; do

		if [ ${txt[$cc]} == ${string}_ENI1.txt ]; then
			stimBeh+="-stim_times $x timing_files/${txt[$cc]} \"BLOCK(0.5,1)\" -stim_label $x beh_${nam[$cc]} "

		elif [ ${txt[$cc]} == ${string}_RI.txt ] || [ ${txt[$cc]} == ${string}_RP.txt ]; then
			stimBeh+="-stim_times $x timing_files/${txt[$cc]} \"BLOCK(6,1)\" -stim_label $x beh_${nam[$cc]} "

		elif [[ ${nam[$cc]} == Int_* ]] || [[ ${nam[$cc]} == Seed* ]]; then
			stimBeh+="-stim_file $x ${txt[$cc]} -stim_label $x ${nam[$cc]} "

		else
			stimBeh+="-stim_times_AM2 $x timing_files/${txt[$cc]} \"dmBLOCK(1)\" -stim_label $x beh_${nam[$cc]} "
		fi

		let x=$[$x+1]
		let cc=$[$cc+1]
	done


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
    -errts ${h_out}_errts \
    -rout -tout" > ${h_out}_deconv.sh
}




### --- Step 0: organize --- ###
#
# determine number of blocks, phases
# determine $input
# set reference variables


cd $ppiDir

# determine input, number of phases/blocks
> tmp.txt
unset input

for a in run*scale+tlrc.HEAD; do

	file=${a%.*}
	input+="$file "

	tmp=${a%_*}
	run=${a%%_*}
	phase=${tmp#*_}

	echo -e "$run \t $phase" >> tmp.txt
done

awk -F '\t' '{print $2}' tmp.txt | sort | uniq -c > phase_list2.txt
rm tmp.txt


# set phase, block arrays
blockArr=(`cat phase_list2.txt | awk '{print $1}'`)
phaseArr=(`cat phase_list2.txt | awk '{print $2}'`)
phaseLen=${#phaseArr[@]}


# set reference variables
firstScale="$(set -- *scale+tlrc.HEAD; echo "$1")"
ref=${firstScale%.*}
TR=`3dinfo -tr $ref`




### --- Step 1: clean data and contrasts --- ###
#
# Create a clean version of each deconvolution by remove effects of no interest.
# Simulate an ideal BOLD response, and generate behavioral contrasts.
#
# The construction of behavior vectors is written for the
# AutismOlfactory dataset, and won't work for other studies


### Make clean data
for i in ${deconList[@]}; do
	if [ ! -f tmp_CleanData_${i}+tlrc.HEAD ]; then


		# Determine effects of no interest - those not in $behInterest
		hold=`grep ColumnLabels X.${i}.xmat.1D | sed 's/\#//g' | sed 's/\;//g'`
		tmp=${hold#*\"}; colLab=(${tmp%\"*})

		unset sel
		for j in ${!colLab[@]}; do
			tmp1=${colLab[$j]#*_}; test=${tmp1%?*}
			MatchString $test "${behInterest[@]}"
			if [ $? == 1 ]; then
				sel+="$j "
			fi
		done


		# remove extra sub-brick added by 3dREMlfit
		sub=`3dinfo -nv ${i}_cbucket_REML+tlrc`
		less=$(($sub-2))
		3dTcat -prefix tmp_${i}_cbucket -tr $TR "${i}_cbucket_REML+tlrc[0..${less}]"


		# create file of unwanted data, replace censored TRs w/mean of neighbors
		3dSynthesize -prefix tmp_effNoInt_AO_${i} -matrix X.${i}.xmat.1D \
		-cbucket tmp_${i}_cbucket+tlrc -select $sel -cenfill nbhr


		# remove unwanted regressors from combined run to produce Clean Data
		3dTcat -prefix tmp_all_runs_${i} -tr $TR $input
		3dcalc -a tmp_all_runs_${i}+tlrc -b tmp_effNoInt_AO_${i}+tlrc -expr 'a-b' -prefix tmp_CleanData_${i}


		# check
		if [ ! -f tmp_CleanData_${i}+tlrc.HEAD ]; then
			echo >&2
			echo "Creating CleanData failed. Exit 1." >&2
			echo >&2
			exit 1
		fi
	fi
done


# Create ideal response function
if [ ! -f ClassicBold.1D ]; then
	waver -dt $TR -GAM -inline 1@1 > ClassicBold.1D
fi


### Create behavior vectors
if [ ! -s Beh_FUBO_contrast.1D ]; then

	cp ${timingDir}/Contrast/${string}_ppi.txt .
	count=`cat ${string}_ppi.txt | wc -l`

	> FBO_contrast.txt
	> UBO_contrast.txt
	> FUBO_contrast.txt


	# Extract ppi.txt info
	for(( a=1; a<=$count; a++ )); do

		value=`sed -n ${a}p ${string}_ppi.txt`


		# FBO & UBO
		if [ $value == 1 ]; then
			echo "1" >> FBO_contrast.txt
			echo "0" >> UBO_contrast.txt
		elif [ $value == -1 ]; then
			echo "0" >> FBO_contrast.txt
			echo "1" >> UBO_contrast.txt
		else
			echo "0" >> FBO_contrast.txt
			echo "0" >> UBO_contrast.txt
		fi


		# FUBO
		if [ $value != 0 ]; then
			echo "1" >> FUBO_contrast.txt
		else
			echo "0" >> FUBO_contrast.txt
		fi
	done

	ConvertDset -o_1D -input FBO_contrast.txt -prefix Beh_FBO_contrast && mv Beh_FBO_contrast.1D.dset Beh_FBO_contrast.1D
	ConvertDset -o_1D -input UBO_contrast.txt -prefix Beh_UBO_contrast && mv Beh_UBO_contrast.1D.dset Beh_UBO_contrast.1D
	ConvertDset -o_1D -input FUBO_contrast.txt -prefix Beh_FUBO_contrast && mv Beh_FUBO_contrast.1D.dset Beh_FUBO_contrast.1D
fi


# CA & Mask
if [ ! -f Beh_Mask_contrast.1D ]; then

	cp ${timingDir}/${string}_TF_* .
	ConvertDset -o_1D -input ${string}_TF_Mask.txt -prefix Beh_Mask_contrast && mv Beh_Mask_contrast.1D.dset Beh_Mask_contrast.1D
	ConvertDset -o_1D -input ${string}_TF_CA.txt -prefix Beh_CA_contrast && mv Beh_CA_contrast.1D.dset Beh_CA_contrast.1D
fi


# Odor
if [ ! -f Beh_Odor_contrast.1D ]; then

	cp ${timingDir}/${string}_CA_Odor.txt .
	cat ${string}_CA_Odor.txt | awk '{print $2}' > tmp_Odor.txt
	ConvertDset -o_1D -input tmp_Odor.txt -prefix Beh_Odor_contrast && mv Beh_Odor_contrast.1D.dset Beh_Odor_contrast.1D
fi


# check
for test in Beh_{FBO,UBO,FUBO,Mask,CA,Odor}_contrast.1D; do
	if [ ! -f $test ]; then
		echo >&2
		echo "Creating $test failed. Exit 2." >&2
		echo >&2
		exit 2
	fi
done




### --- Step 2: extract seed time series --- ###
#
# Construct seeds, and extract mean time series (TS) from Clean Data.
# TS is deconvolved for BOLD response, rendering "neural" TS.
#
# Functional TR = 2000 ms, stimulus duration is measured in 100 ms,
# so the TS is upsampled in order to extract "neural" TS associated
# w/the behavior.
#
# Then the TS is downsampled and reconvolved with BOLD to generate
# a task-associated BOLD TS.
#
# Much of this is specific to the AutismOlfactory study


# Get seed TS from e/deconvolution
c=0; while [ $c -lt $seedLen ]; do


	# Make seeds
	seed=Seed_${seedName[$c]}

	if [ ! -f ${seed}+tlrc.HEAD ]; then
		echo ${seedCoord[$c]} > tmp_${seed}.txt
		3dUndump -prefix $seed -master $ref -srad 3 -xyz tmp_${seed}.txt
	fi

	for j in ${deconList[@]}; do
		if [ ! -f tmp_HRes_${seed}_${j}_neural.1D ] && [ ! -f ${seed}_${j}_TS_CA.1D ]; then


			# Get seed TS from each CleanData, solve RHS for neural
			3dmaskave -quiet -mask ${seed}+tlrc tmp_CleanData_${j}+tlrc > ${seed}_${j}_timeSeries.1D
			3dTfitter -RHS ${seed}_${j}_timeSeries.1D -FALTUNG ClassicBold.1D tmp_${seed}_${j}_neural 012 0


			# Resample seed TS from 2s to 0.1s resolution (HRes)
			> tmp_HRes_${seed}_${j}.txt
			1dtranspose tmp_${seed}_${j}_neural.1D > tmp_Trans_${seed}_${j}_neural.1D
			data=(`cat tmp_Trans_${seed}_${j}_neural.1D`)
			dataLen=${#data[@]}

			cc=0; while [ $cc -lt $dataLen ]; do
				for((i=1; i<=20; i++)); do
					echo ${data[$cc]} >> tmp_HRes_${seed}_${j}.txt
				done
				let cc=$[$cc+1]
			done


			# Create 1D file
			ConvertDset -o_1D -input tmp_HRes_${seed}_${j}.txt -prefix tmp_HRes_${seed}_${j}_neural
			mv tmp_HRes_${seed}_${j}_neural.1D.dset tmp_HRes_${seed}_${j}_neural.1D
		fi
	done

	let c=$[$c+1]
done


### Extract behavior TS from appropriate deconv from e/seed
for i in ${deconList[@]}; do


	# pull nested array info
	conList=($(eval echo \${arr${i}[@]}))

	for j in ${conList[@]}; do
		for k in ${seedName[@]}; do
			if [ ! -f Seed_${k}_${i}_TS_${j}.1D ]; then


				# Extract seed beh neural timeseries, resample back to 2s (LRes)
				1deval -a tmp_HRes_Seed_${k}_${i}_neural.1D -b Beh_${j}_contrast.1D -expr 'a*b' > tmp_Seed_${k}_${i}_neural_${j}_beh.1D
				cat tmp_Seed_${k}_${i}_neural_${j}_beh.1D | awk -v n=20 'NR%n==0' > tmp_LRes_Seed_${k}_${i}_neural_${j}.txt
				ConvertDset -o_1D -input tmp_LRes_Seed_${k}_${i}_neural_${j}.txt -prefix tmp_LRes_Seed_${k}_${i}_neural_${j}
				mv tmp_LRes_Seed_${k}_${i}_neural_${j}.1D.dset tmp_LRes_Seed_${k}_${i}_neural_${j}.1D


				# add BOLD back to TS
				num=`cat tmp_Trans_Seed_${k}_${i}_neural.1D | wc -l`
				waver -GAM -peak 1 -dt $TR -input tmp_LRes_Seed_${k}_${i}_neural_${j}.1D -numout $num > Seed_${k}_${i}_TS_${j}.1D


				# check
				if [ ! -f Seed_${k}_${i}_TS_${j}.1D ]; then
					echo >&2
					echo "Creating Seed_${k}_${i}_TS_${j}.1D failed. Exit 3." >&2
					echo >&2
					exit 3
				fi
			fi
		done
	done
done




### --- Step 3: Deconvolve --- ###
#
# A deconvolution script (foo_deconv.sh) is generated and ran for
# each planned deconvolution.
#
# This is a stripped version of our normal GenDecon call syntax, and
# has been adjusted for the AutismOlfactory PPI.
#
# REML will be run to extract correlation matrix.


# for each conducted decon & seed
c=0; while [ $c -lt $phaseLen ]; do
	for i in ${deconList[@]}; do
		for j in ${seedName[@]}; do


			# extract nested array info
			holdName=($(eval echo \${nam${i}[@]}))
			holdTxt=($(eval echo \${txt${i}[@]}))


			# append holdName/Txt to holds Seed_TS.1D info
			holdName+=(Seed_${j})
			holdTxt+=(Seed_${j}_${i}_timeSeries.1D)

			conList=($(eval echo \${arr${i}[@]}))

			for k in ${conList[@]}; do
				if [ -f Seed_${j}_${i}_TS_${k}.1D ]; then

					holdTxt+=(Seed_${j}_${i}_TS_${k}.1D)
					holdName+=(Int_$k)
				fi
			done


			if [ ! -f X.PPI_${i}_${j}.xmat.1D ]; then


				# generate script
				GenDecon ${phaseArr[$c]} ${blockArr[$c]} "$input" PPI_${i}_$j ${#holdName[@]} ${holdName[@]} ${holdTxt[@]}


				# run script
				if [ -f PPI_${i}_${j}_stats.REML_cmd ]; then
					rm PPI_${i}_${j}_stats.REML_cmd
				fi
				source PPI_${i}_${j}_deconv.sh


				# check
				if [ ! -f X.PPI_${i}_${j}.xmat.1D ]; then
					echo >&2
					echo "Creating Decon failed for ${i}_${j}. Exit 4." >&2
					echo >&2
					exit 4
				fi
			fi


			# run REML script
			if [ ! -f PPI_${i}_${j}_stats_REML+tlrc.HEAD ]; then
				tcsh -x PPI_${i}_${j}_stats.REML_cmd -dsort ${phaseArr[$c]}_WMe_rall+tlrc
			fi


			# kill if REMl failed
			if [ ! -f PPI_${i}_${j}_stats_REML+tlrc.HEAD ]; then
				echo "" >&2
				echo "REML failed on ${i}_${j}. Exit 5." >&2
				echo "" >&2
				exit 5
			fi
		done
	done
	let c=$[$c+1]
done
