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
    D)  DRYRUN='1';;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

# list scripts to run
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
SCRIPTS=("${SCRIPT_DIR}/exp"*)

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
        --mem=6GB
        --account='daniel_marcus'
        --partition='tier2_cpu'
        --exclude=highmem01,highmem02
        call_python.sh $script $SHORT)
    echo ""
    echo "${COMMAND[@]}"

    if [[ $DRYRUN == '0' ]]
    then
        "${COMMAND[@]}"
    fi

done
