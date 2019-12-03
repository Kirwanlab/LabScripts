#!/bin/bash

#SBATCH --time=15:00:00   # walltime
#SBATCH --ntasks=6   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=8gb   # memory per CPU core
#SBATCH -J "RS2"   # job name

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
# This script does single-subject regression for each participant.
# 	A correlation matrix is then projected out using the updated anaticor
#	method (see example 11 of afni_proc.py, and https://www.sciencedirect.com/science/article/pii/S1053811910006580).
#
# Desired output file is errts_fanaticor+tlrc.
#
# Also models noise (ACF) for potential Monte Carlo simulations, and
#	generates a review file (out_summary.txt). This review file is
#	used in subsequent steps.
#
# Major steps following the comment ### --- Foo --- ###,
# 	minor steps are also annotated.


subj=$1
phase=$2


### --- Set Up --- ###

multiBand=1				# for vibration patch

workDir=~/compute/MNI_NREM/MRI_processed/derivatives/${subj}/ses-$phase




### --- Motion Files --- ###
#
# Construct de-meaned and derivative motion files.
# Also create a censor file (motion*outlier).
#
# Mutliband patch - a vibration may be occuring from
# the protocol, resulting in a large number of volumes
# being excluded. Volreg was used to account for this
# vibration. Currently trying the decon with only the
# outlier censoring.


cd $workDir

if [ ! -s censor_Resting_${phase}_combined_2.1D ]; then

	cat dfile_Resting_${phase}.1D > dfile_rall.1D
	cat out.cen.RS.1D > outcount_Resting_${phase}_censor.1D

	# demean
	1d_tool.py -infile dfile_rall.1D -set_nruns 1 -demean -write motion_demean.1D
	1d_tool.py -infile motion_demean.1D -set_nruns 1 -split_into_pad_runs mot_demean

	# derivative
	1d_tool.py -infile dfile_rall.1D -set_nruns 1 -derivative -demean -write motion_deriv.1D
	1d_tool.py -infile motion_deriv.1D -set_nruns 1  -split_into_pad_runs mot_deriv

	# censor - 0.2 <- be more conservative (than 0.3) for RS analyses
	1d_tool.py -infile dfile_rall.1D -set_nruns 1 -show_censor_count -censor_prev_TR -censor_motion 0.2 motion_Resting_${phase}

	if [ $multiBand == 1 ]; then
		cp outcount_Resting_${phase}_censor.1D censor_Resting_${phase}_combined_2.1D
	else
		1deval -a motion_Resting_${phase}_censor.1D -b outcount_Resting_${phase}_censor.1D -expr "a*b" > censor_Resting_${phase}_combined_2.1D
	fi
fi


# note TRs that were not censored
ktrs=`1d_tool.py -infile censor_Resting_${phase}_combined_2.1D -show_trs_uncensored encoded`




### --- Deconvolve --- ###
#
# Run a PCA on cleaned data to remove CSF signal in order to
# collapse that signal into baseline. Produce residual time
# series matrix that has been cleaned by anaticor model

tr_counts=`3dinfo -ntimes Resting_${phase}+orig`
len_tr=`3dinfo -tr Resting_${phase}+orig`
pol_time=$(echo $(echo $tr_counts*$len_tr | bc)/150 | bc -l)
pol=$((1 + `printf "%.0f" $pol_time`))

if [ ! -s ROIPC_Svent.1D ]; then

	# set up, run PCA
	1d_tool.py -set_run_lengths $tr_counts -select_runs 1 -infile censor_Resting_${phase}_combined_2.1D -write tmp_censor.1D
	3dTproject -polort $pol -prefix tmp_det_pcin -censor tmp_censor.1D -cenmode KILL -input tmp_epi_clean+tlrc
	3dpc -mask final_mask_CSF_eroded+tlrc -pcsave 3 -prefix tmp_ROIPC_Svent tmp_det_pcin+tlrc

	# write PCA files for 3dDeconvolve
	1d_tool.py \
	-censor_fill_parent tmp_censor.1D \
	-infile tmp_ROIPC_Svent_vec.1D \
	-write - | 1d_tool.py \
	-set_run_lengths $tr_counts \
	-pad_into_many_runs 1 1 \
	-infile - -write ROIPC_Svent.1D
fi


# deconvolve
if [ ! -f X.xmat.1D ]; then

	3dDeconvolve -input Resting_${phase}_scale+tlrc \
	-censor censor_Resting_${phase}_combined_2.1D \
	-ortvec ROIPC_Svent.1D ROIPC_Svent_Resting_${phase} \
	-polort A \
	-num_stimts 12 \
	-stim_file 1 mot_demean.r01.1D'[0]' -stim_base 1 -stim_label 1 roll_01  \
	-stim_file 2 mot_demean.r01.1D'[1]' -stim_base 2 -stim_label 2 pitch_01 \
	-stim_file 3 mot_demean.r01.1D'[2]' -stim_base 3 -stim_label 3 yaw_01   \
	-stim_file 4 mot_demean.r01.1D'[3]' -stim_base 4 -stim_label 4 dS_01    \
	-stim_file 5 mot_demean.r01.1D'[4]' -stim_base 5 -stim_label 5 dL_01    \
	-stim_file 6 mot_demean.r01.1D'[5]' -stim_base 6 -stim_label 6 dP_01    \
	-stim_file 7 mot_deriv.r01.1D'[0]' -stim_base 7 -stim_label 7 roll_02   \
	-stim_file 8 mot_deriv.r01.1D'[1]' -stim_base 8 -stim_label 8 pitch_02  \
	-stim_file 9 mot_deriv.r01.1D'[2]' -stim_base 9 -stim_label 9 yaw_02    \
	-stim_file 10 mot_deriv.r01.1D'[3]' -stim_base 10 -stim_label 10 dS_02  \
	-stim_file 11 mot_deriv.r01.1D'[4]' -stim_base 11 -stim_label 11 dL_02  \
	-stim_file 12 mot_deriv.r01.1D'[5]' -stim_base 12 -stim_label 12 dP_02  \
    -fout -tout -x1D X.xmat.1D -xjpeg X.jpg                                 \
    -x1D_uncensored X.nocensor.xmat.1D                                      \
    -fitts fitts   		                                                    \
    -errts errts    		                                                \
    -x1D_stop                                                               \
    -bucket stats


	if [ ! -f X.xmat.1D ]; then
		echo >&2
		echo "Deconvolve failed - no X.xmat.1d detected. Exit 1" >&2
		echo >&2; exit 1
	fi

	# display any large pairwise correlations from the X-matrix
	1d_tool.py -show_cormat_warnings -infile X.xmat.1D | tee out_cor_matrix_warn.txt

	# baseline regressors
	reg_cols=`1d_tool.py -infile X.nocensor.xmat.1D -show_indices_interest`
	3dTstat -sum -prefix sum_ideal.1D X.nocensor.xmat.1D"[$reg_cols]"
	1dcat X.nocensor.xmat.1D"[$reg_cols]" > X.stim.xmat.1D
fi


# Project out regression matrix
if [ ! -f errts_tproject+tlrc.HEAD ];then

	3dTproject \
	-polort 0 \
	-input Resting_${phase}_scale+tlrc.HEAD \
	-censor censor_Resting_${phase}_combined_2.1D \
	-cenmode ZERO \
	-ort X.nocensor.xmat.1D \
	-prefix errts_tproject
fi


# Anaticor
if [ ! -f errts_fanaticor+tlrc.HEAD ]; then

	3dTcat -prefix tmp_all_Resting_${phase} Resting_${phase}_scale+tlrc
	3dTcat -prefix tmp_all_volreg tmp_epi_clean+tlrc
	3dcalc -a tmp_all_volreg+tlrc -b final_mask_CSF_eroded+tlrc \
	-expr "a*bool(b)" -datum float -prefix tmp_all_volreg_mask

	3dmerge -1blur_fwhm 30 -doall -prefix Local_WMe_rall tmp_all_volreg_mask+tlrc

	3dTproject \
	-polort 0 \
	-input Resting_${phase}_scale+tlrc \
	-censor censor_Resting_${phase}_combined_2.1D \
	-cenmode ZERO \
	-dsort Local_WMe_rall+tlrc \
	-ort X.nocensor.xmat.1D \
	-prefix errts_fanaticor

	if [ ! -f errts_fanaticor+tlrc.HEAD ]; then
		echo >&2
		echo "Anaticor step failed. Exit 2" >&2
		echo >&2; exit 2
	fi
fi




### --- Model Noise, Generate Review --- ###

if [ ! -f TSNR+tlrc.HEAD ]; then

	# SNR
	3dTstat -mean -prefix tmp_signal_all tmp_all_Resting_${phase}+tlrc"[$ktrs]"
	3dTstat -stdev -prefix tmp_noise_all errts_fanaticor+tlrc"[$ktrs]"

	3dcalc \
	-a tmp_signal_all+tlrc \
	-b tmp_noise_all+tlrc \
	-c mask_epi_anat+tlrc \
	-expr 'c*a/b' \
	-prefix TSNR


	#GCOR
	3dTnorm -norm2 -prefix tmp_errts_unit errts_fanaticor+tlrc
	3dmaskave -quiet -mask full_mask+tlrc tmp_errts_unit+tlrc > gmean_errts_unit.1D
	3dTstat -sos -prefix - gmean_errts_unit.1D\' > out_gcor.1D


	# compute corr
	3dcalc -a tmp_errts_unit+tlrc -b gmean_errts_unit.1D -expr 'a*b' -prefix tmp_fm_WMe
	3dmaskave -q -mask tmp_fm_WMe+tlrc tmp_errts_unit+tlrc > mean_unit_WMe.1D
	3dcalc -a tmp_errts_unit+tlrc -b mean_unit_WMe.1D -expr 'a*b' -prefix tmp_DP_WMe
	3dTstat -sum -prefix corr_af_WMe tmp_DP_WMe+tlrc

	3dcalc -a final_mask_CSF_eroded+tlrc -b full_mask+tlrc -expr 'a*b' -prefix tmp_fm_Svent
	3dmaskave -q -mask tmp_fm_Svent+tlrc tmp_errts_unit+tlrc > mean_unit_Svent
	3dcalc -a tmp_errts_unit+tlrc -b mean_unit_Svent.1D -expr 'a*b' -prefix tmp_DP_Svent
	3dTstat -sum -prefix corr_af_Svent tmp_DP_Svent_tlrc
fi


# noise estimation
if [ ! -s blur_errts.1D ]; then

	> blur_est.1D
	> blur_epits.1D
	> blur_errts.1D

	mkdir files_ACF
	trs=`1d_tool.py -infile X.xmat.1D -show_trs_uncensored encoded -show_trs_run 1`

	3dFWHMx \
	-mask mask_epi_anat+tlrc \
	-ACF files_ACF/out_3dFWHMx_ACF_epits.1D \
	tmp_all_Resting_${phase}+tlrc"[$trs]" >> blur_epits.1D

	3dFWHMx \
	-mask mask_epi_anat+tlrc \
	-ACF files_ACF/out_3dFWHMx_ACF_errts.1D \
	errts_fanaticor+tlrc"[$trs]" >> blur_errts.1D

	if [ ! -s blur_errts.1D ]; then
		echo >&2
		echo "3dFWHMx failed. Exit 3" >&2
		echo >&2; exit 3
	fi

	# averge FWHM blur
	blurE0=`3dTstat -mean -prefix - blur_epits.1D'{0..$(2)}'\'`
	blurE1=`3dTstat -mean -prefix - blur_epits.1D'{1..$(2)}'\'`
	echo "$blurE0     #epits FWHM blur est" >> blur_est.1D
	echo "$blurE1     #epits ACF blur est" >> blur_est.1D

	blurR0=`3dTstat -mean -prefix - blur_errts.1D'{0..$(2)}'\'`
	blurR1=`3dTstat -mean -prefix - blur_errts.1D'{1..$(2)}'\'`
	echo "$blurR0     #errts FWHM blur est" >> blur_est.1D
	echo "$blurR1     #errts ACF blur est" >> blur_est.1D
fi


# Review scripts
if [ ! -s out_summary.txt ]; then

	3dcopy Resting_${phase}+orig pb00.${subj}.r01.tcat
	3dcopy tmp_epi_clean+tlrc pb02.${subj}.r01.volreg

	gen_ss_review_scripts.py \
	-subj ${subj} \
	-rm_trs 0 \
	-motion_dset dfile_rall.1D \
	-outlier_dset outcount.RS.1D \
	-enorm_dset  motion_Resting_${phase}_enorm.1D \
	-mot_limit 0.2 \
	-out_limit 0.05 \
	-xmat_regress X.xmat.1D \
	-xmat_uncensored X.nocensor.xmat.1D \
	-errts_dset errts_fanaticor+tlrc \
	-final_anat final_anat+tlrc \
	-final_view tlrc \
	-exit0

	./\@ss_review_basic | tee out_summary.txt
fi
