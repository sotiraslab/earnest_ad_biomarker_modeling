# list scripts to run
SCRIPTS=('plt0_boxplot_all_individual_predictors_global_cog.py'
         'plt1_boxlot_combo_models_vs_covariates.py'
         'plt2_boxplot_combo_models_vs_binary.py'
         'plt3_combine_binary_baseline_experiments_neuropsych.py'
)

echo ""
echo "RUNNING PLOTTING SCRIPTS"
echo "------------------------"

for script in "${SCRIPTS[@]}"
do  
    COMMAND=(python ${script})
    echo ""
    echo "==================================================" 
    echo "* Beginning: ${script}"
    echo "* Command: ${COMMAND[@]}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    "${COMMAND[@]}"
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "* Done!"
    echo "=================================================="
done