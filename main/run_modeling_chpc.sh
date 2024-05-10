#!//bin/bash

# Define help text for the script.
usage() {
	echo ""
	echo "Usage:	run_modeling_chpc.sh [-S -D]"
	echo ""
	echo "Description:"
	echo "    Run all cross-validated modeling experiments by submitting SLURM jobs."
	echo "Arguments:"
    echo "    -S Short.  Limits hyperparameter tuning for all SVM experiments.  No difference for linear models."
	echo "    -D Dryrun.  Commands are printed, but jobs are not submitted"
    echo ""
	}

# read arguments
SHORT=''
DRYRUN='0'
while getopts hSD arg
do
	case $arg in
	h)	usage
		exit 0;;
	S)	SHORT='--short';;
    D)  DRYRUN=='1';;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

# list scripts to run
SCRIPTS=('exp0_individual_atn_models_global_cognition.py'
         'exp1_svms_global_cognition.py'
         'exp2_combo_atn_models_global_cognition.py'
         'exp3_combo_atn_models_global_cognition_vs_binary.py'
         'exp4a_combo_atn_models_memory_vs_binary.py'
         'exp4b_svms_memory.py'
         'exp5a_combo_atn_models_executive_functioning_vs_binary.py'
         'exp5b_svms_executive_functioning.py'
         'exp6a_combo_atn_models_language_vs_binary.py'
         'exp6b_svms_language.py'
         'exp7a_combo_atn_models_visuospatial_vs_binary.py'
         'exp7b_svms_visuospatial.py'
         'exp8a_combo_atn_models_longitudinal_cognition_vs_binary.py'
         'exp8b_svms_longitudinal_cognition.py'
         'exp9a_preclinical_combo_atn_models.py'
         'exp9b_preclinical_combo_atn_models_vs_binary.py'
         'exp9c_preclinical_svms.py'
         'exp10a_test_all_predictors_with_csf.py'
         'exp10b_test_all_predictors_with_csf_vs_baseline.py'
)

# silently load python on CHPC
module load python
source ~/.bashrc
source ~/miniconda/bin/activate atn_cognition

for script in "${SCRIPTS[@]}"
do
    name=$(basename $script .py)
    mkdir -p "logs"
    log="logs/${name}.log"
    
    if [[ $SHORT == '--short' ]]
    then
        log="logs/${name}_short.log"
    fi

    COMMAND=(sbatch
        -J $name
        -o $log
        -t '72:00:00'
        -N 1
        -n 1
        --mem=128GB
        --account='aristeidis_sotiras'
        --partition='tier2_cpu'
        --exclude=highmem01,highmem02
        call_python.sh $script $SHORT)
    echo "${COMMAND[@]}"

    if [[ $DRYRUN == '0' ]]
    then
        "${COMMAND[@]}"
    fi

done
