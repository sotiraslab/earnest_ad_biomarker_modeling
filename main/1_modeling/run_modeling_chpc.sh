#!//bin/bash

# silently load python on CHPC
module load python
source activate atn_cognition

SCRIPTS=('exp0_individual_atn_models_global_cognition.py'
         'exp1_svms_global_cognition.py'
)

for script in "${SCRIPTS[@]}"
do
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
        call_python.sh $script)
    echo "${COMMAND[@]}"
    "${COMMAND[@]}"

done
