#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import matplotlib.pyplot as plt
import pandas as pd

from common import load_results, set_labels_baseline_exp
from atn_modeling.helpers import results_boxplot

# load data
lm_CU = load_results('exp9a_CU_combo_atn_models', 'results.csv')
svm_CU = load_results('exp9c_CU_svms', 'results.csv')
concat_CU = pd.concat([lm_CU, svm_CU])
concat_CU = concat_CU.loc[~concat_CU['model'].str.contains('PVC')]

lm_CI = load_results('exp9d_CI_combo_atn_models', 'results.csv')
svm_CI = load_results('exp9f_CI_svms', 'results.csv')
concat_CI = pd.concat([lm_CI, svm_CI])
concat_CI = concat_CI.loc[~concat_CI['model'].str.contains('PVC')]

# general resources
palette = (['gray'] +
    ['#FFDDAA'] * 3 +
    ['#F7A934'] +
    ['#E899EE'] * 3 +
    ['#B74CBF'] +
    ['#F29D9D'] * 3 +
    ['#FC4646'] +
    ['#A5DBF2'] * 3 +
    ['#08A3E5'])


# plot CU
plot_path = os.path.join('figures', 'boxplot_combo_models_vs_covariates_CU.svg')
n_train = concat_CU['ntrain'].values[0]
n_test = concat_CU['ntest'].values[0]
fig, stats_CU = results_boxplot(concat_CU, groupby='model', baseline='Baseline',
                                stats_vs_baseline=True, palette=palette,
                                n_train=n_train, n_test=n_test, font_file='arial.ttf')
set_labels_baseline_exp(fig)

plt.title('Clinically Unimpaired')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)

# plot CI
plot_path = os.path.join('figures', 'boxplot_combo_models_vs_covariates_CI.svg')
n_train = concat_CI['ntrain'].values[0]
n_test = concat_CI['ntest'].values[0]
fig, stats_CI = results_boxplot(concat_CI, groupby='model', baseline='Baseline',
                                stats_vs_baseline=True, palette=palette,
                                n_train=n_train, n_test=n_test, font_file='arial.ttf')
set_labels_baseline_exp(fig)

plt.title('Clinically Impaired')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)