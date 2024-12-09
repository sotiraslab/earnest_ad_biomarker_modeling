# Evaluating biomarker variability in Alzheimer's Disease
This repository contain's code for research from the [Sotiras Lab](https://www.mir.wustl.edu/employees/aristeidis-sotiras/) investigating how different imaging biomarker operationalizations affect cognitive prediction in AD.

## Citation

This work is currently shared in the following preprint: **Comprehensive evaluation of AT(N) imaging biomarkers for predicting cognition** ([http://dx.doi.org/10.1101/2024.11.25.24317943](http://dx.doi.org/10.1101/2024.11.25.24317943)).

## User Guide

This section will document how this repository is organized and how to use the code to reproduce analyses.

### 1. Acquire input data (ADNI)

You must be approved to access ADNI data ([see here](https://adni.loni.usc.edu/data-samples/adni-data/)).  From the [ADNI data repository website](https://ida.loni.usc.edu/login.jsp), you need to download the following tables:

- UC Berkeley amyloid 6mm analysis (`UCBERKELEY_TAU_6MM_[date].csv`)
- UC Berkeley tau 6mm analysis (`UCBERKELEY_TAU_6MM_[date].csv`)
- UC Berkeley tau 6mm analysis with partial volume correction(`UCBERKELEY_TAUPVC_6MM_[date].csv`)
- ADSP Phenotype Harmonization Scores (`ADSP_PHC_COGN_[date].csv`)

**All these should be placed in `earnest_ad_biomarker_modeling/main/inputs`.**  They will only be searched for in this folder.

You also need to download and install the ADNIMERGE R package, whose source code is also available from the ADNI data repository.  There are some instructions for [how to install this here](https://adni.bitbucket.io/).

### 2. Install required software

There are a few pieces of software that need to be installed:

- **R packages**: In addition to ADNIMERGE R (see previous section), you need to have the following R packages installed:

  - colormap
  - ggplot2
  - ggseg
  - gsubfn
  - lme4
  - lubridate
  - stringr
  - svglite
  - tableone
  - this.path
  - tidyverse

- **Python code**: This repository contains a Python package called atn_modeling.  Installing it with pip should install all the required Python dependencies (matplotlib, numpy, pandas, pingouin, scipy, scikit-learn>=1.4, statsmodels).  Here is how it can be installed:

  ```
  git clone https://github.com/sotiraslab/earnest_ad_biomarker_modeling
  cd earnest_ad_biomarker_modeling
  pip install .
  ```

### 3. Run analyses

All of the code to generate analyses are in the `main` folder.  Directly in the folder are several R and Python scripts which do the main work:

1. `build_datasets.R` will generate all the tables of subjects which are used for modeling.  These are placed in `main/datasets`.
2. All the scripts starting `exp[...].py` are "experiments", which are cross-validated modeling analyses, both with linear regressions and SVMs.  Parameters from these models will be saved in the `main/outputs` folder.  All these files can be run as scripts, e.g. `python exp0_individual_atn_models_global_cognition.py` from the command line.
   - **NOTE**: SVM experiments can take a long time to run, with hyperparmeter grid search (several hours).  This is required in order to full replicate the analyses presented in the paper.  For quicker results, all experiments scripts accept a `--short` flag, which will limit the hyperparameter search space.  E.g., `python experiment1_svms_global_cognition.py --short`.  
3. All the scripts starting `plt[...].py` and `plt[...].R` are scripts which generate figures derived from experimental results.  They will save outputs in the `main/figures` folder.

Note that these scripts **need to generally be run in the order just listed**.  That is, first run `build_datasets.R` to generate the tables of data used for analysis.  Then, run any cross-validation experiment scripts.  Finally, generate figures with any of the plotting scripts.  However, not all plots need all experiments to be run in order to be generated.

Provided you have run all experiments, you can use `run_all_plots.sh` to generate all main text & supplemental figures.  There is not a script to run all experiments iteratively, since doing so serially would be prohibitively long.  But there are scripts for running on a computing cluster in parallel (which could be edited to run things serially).

#### Cluster computing

The authors ran these analyses using the [Center for High Performance Computing cluster at Washington University in St Louis](https://www.mir.wustl.edu/research/core-resources/research-computing-and-informatics-facility/).  This allowed for parallel processing of different cross-validation experiments, greatly lowering the time to complete a round of analysis.  There are a couple BASH scripts that demonstrate how this was implemented.

- `main/run_modeling_chpc.sh` is what I used to run cross-validation experiments on a SLURM cluster, each submitted as a separate job.  This scripts also accept additional arguments to have a dry run (`-D`) and to run with quicker hyperparameter search (`-S`).
- `main/run_selection_modeling_chpc.sh` is basically the same, but allows you to select a subset of experiments to run.
- `upload_chpc.py`: Upload the datasets/inputs folder to the cluster.
- `download_chpc.py`: Download results from modeling run on the clustering.

Note that some paths and SLURM options are specific to the author, and will need to be updated.  Please get in contact if you need help amending the code.

## Troubleshooting

Please raise an [issue](https://github.com/sotiraslab/earnest_ad_biomarker_modeling/issues) if you run into an problems or questions using the code in this repository.
