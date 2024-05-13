#!//bin/bash

# Define help text for the script.
usage() {
	echo ""
	echo "Usage:	run_modeling_chpc.sh -i "script.py" [-i "another_script.py"][-S -D]"
	echo ""
	echo "Description:"
	echo "    Run a manually input selection of cross-validated modeling experiments by submitting SLURM jobs."
	echo "Arguments:"
    echo "    -i Input.  A single script to run.  Multiple can be entered as separate argments."
    echo "    -S Short.  Limits hyperparameter tuning for all SVM experiments.  No difference for linear models."
	echo "    -D Dryrun.  Commands are printed, but jobs are not submitted"
    echo ""
	}

# read arguments
SHORT=''
DRYRUN='0'
SCRIPTS=()
while getopts hSDi: arg
do
	case $arg in
	h)	usage
		exit 0;;
	S)	SHORT='--short';;
    D)  DRYRUN='1';;
    i)  SCRIPTS+=("$OPTARG");;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

# check something provide
if [ ${#SCRIPTS[@]} -eq 0 ]
then
    echo ""
    echo "ERROR: No scripts provided!  Exiting"
    usage
    exit 1
fi

# list scripts to run
echo ""
echo "You have entered the following scripts to run:"
echo ""
for val in "${SCRIPTS[@]}"; do
    echo " - $val"
done

echo ""
read -p "Continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

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
