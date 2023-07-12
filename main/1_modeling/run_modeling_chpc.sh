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
)

for script in "${SCRIPTS[@]}"
do
    name=$(basename $script .py)
    log=${name}.slurmlog

    COMMAND=(sbatch -J $name -o $log -t '24:00:00' -N 1 -n 1 --mem=64GB call_python.sh $script)
    echo "${COMMAND[@]}"
    "${COMMAND[@]}"

done
