#!//bin/bash

# silently load python on CHPC
module load python
source activate atn_cognition

SCRIPTS=('control_baseline.py'
'control_baseline_CN.py'
'control_baseline_DEM.py'
'binary_baseline.py'
'binary_baseline_CN.py'
'binary_baseline_DEM.py'
'control_baseline_longitudinal.py'
'binary_baseline_longitudinal.py'
'data_driven_v_established.py'
)

for script in "${SCRIPTS[@]}"
do
    name=$(basename $script .py)
    log=${name}.slurmlog

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
