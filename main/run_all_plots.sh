# list scripts to run
SCRIPTS=('plt0_boxplot_all_individual_predictors_global_cog.py'
         'plt1_boxlot_combo_models_vs_covariates.py'
         'plt2_boxplot_combo_models_vs_binary.py'
         'plt3_combine_binary_baseline_experiments_neuropsych.py'
         'plt4_boxlot_combo_models_vs_covariates_CU_CI.py'
         'plt5_boxplot_combo_models_vs_binary_CU_CI.py'
         'plt6_cutoff_distribution.py'
         'plt7_linear_model_coefficients_distribution.py'
         'plt8_feature_selection_piecharts.py'
         'plt9_taupvc_individual_models.py'
         'plt10_taupvc_combination_models.py'
         'plt12_boxplot_all_individual_predictors_with_csf.py'
         'plt13_boxplot_combination_models_with_csf.py'
         'plt14_boxlot_combo_vs_covariates_longitudinal.py'
         'plt15_boxlot_combo_vs_binary_longitudinal.py'
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

# run Haufe viz scripts separately, as they involve R
script='plt11a_extract_haufe_weights.py'
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

script='plt11b_plot_haufe_weights.R'
COMMAND=(Rscript ${script})
echo ""
echo "==================================================" 
echo "* Beginning: ${script}"
echo "* Command: ${COMMAND[@]}"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
"${COMMAND[@]}"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "* Done!"
echo "=================================================="