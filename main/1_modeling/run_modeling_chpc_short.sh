#!//bin/bash

# silently load python on CHPC
module load python
source activate atn_cognition

SCRIPTS=('exp0_individual_atn_models_global_cognition.py'
         'exp1_svms_global_cognition.py'
)

FLAGS=(''
       '--short')

# use for loop to read all values and indexes
for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
    script=${SCRIPTS[$i]}
    flags=${FLAGS[$i]}

    name=$(basename $script .py)
    log="../../outputs/${name}/slurm.log"

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
        call_python.sh $script $flags)
    echo "${COMMAND[@]}"
    "${COMMAND[@]}"

done