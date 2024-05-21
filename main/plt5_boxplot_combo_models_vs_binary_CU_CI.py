#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import matplotlib.pyplot as plt
import pandas as pd

from atn_modeling.helpers import results_boxplot
from common import load_results

# load data
lm_CU = load_results('exp9b_CU_combo_atn_models_vs_binary', 'results.csv')
svm_CU = load_results('exp9c_CU_svms', 'results.csv')
concat_CU = pd.concat([lm_CU, svm_CU])
concat_CU = concat_CU.loc[~concat_CU['model'].str.contains('PVC')]

lm_CI = load_results('exp9e_CI_combo_atn_models_vs_binary', 'results.csv')
svm_CI = load_results('exp9f_CI_svms', 'results.csv')
concat_CI = pd.concat([lm_CI, svm_CI])
concat_CI = concat_CI.loc[~concat_CI['model'].str.contains('PVC')]

# output
plot_path = os.path.join('figures', os.path.splitext(os.path.basename(__file__))[0] + '.svg')

# general resources
palette = (['#F7A934'] +
    ['#E899EE'] * 3 +
    ['#B74CBF'] +
    ['#F29D9D'] * 3 +
    ['#FC4646'] +
    ['#A5DBF2'] * 3 +
    ['#08A3E5'])

# plot CU
plot_path = os.path.join('figures', 'boxplot_combo_models_vs_binary_CU.svg')
n_train = concat_CU['ntrain'].values[0]
n_test = concat_CU['ntest'].values[0]
fig, _ = results_boxplot(concat_CU, groupby='model', baseline='All binary',
                         stats_vs_baseline=True, palette=palette,
                         n_train=n_train, n_test=n_test, font_file='arial.ttf')

plt.title('Clinically Unimpaired')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)

# plot CI
plot_path = os.path.join('figures', 'boxplot_combo_models_vs_binary_CI.svg')
n_train = concat_CI['ntrain'].values[0]
n_test = concat_CI['ntest'].values[0]
fig, _ = results_boxplot(concat_CI, groupby='model', baseline='All binary',
                         stats_vs_baseline=True, palette=palette,
                         n_train=n_train, n_test=n_test, font_file='arial.ttf')

plt.title('Clinically Impaired')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)
