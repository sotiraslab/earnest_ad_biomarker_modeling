#!//bin/bash

# silently load python on CHPC
module load python
source activate atn_cognition

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
)


FLAGS=(''
       '--short'
       ''
       ''
       ''
       '--short'
       ''
       '--short'
       ''
       '--short'
       ''
       '--short'
       ''
       '--short')

# use for loop to read all values and indexes
for (( i=0; i<${#SCRIPTS[@]}; i++ ));
do
    script=${SCRIPTS[$i]}
    flags=${FLAGS[$i]}

    name=$(basename $script .py)
    mkdir -p "../../outputs/logs"
    log="../../outputs/logs/${name}_short.log"

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
        call_python.sh $script --rerun $flags)
    echo "${COMMAND[@]}"
    "${COMMAND[@]}"

done