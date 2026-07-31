#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os

import matplotlib.pyplot as plt
import pandas as pd

from atn_modeling.helpers import results_boxplot
from common import load_results, set_labels_baseline_exp, set_labels_binary_exp

#%% Plot combination vs baseline

exp_lm = load_results('exp2_comboVSsolo_MMSE', 'results.csv')
exp_svm = load_results('exp1_svm_MMSE', 'results.csv')
spearman_lm = pd.read_csv(os.path.join('figures', 'spearman_correlations', 'combo_vs_baseline.csv'))
spearman_svm = pd.read_csv(os.path.join('figures', 'spearman_correlations', 'svm.csv'))

concat = pd.concat([exp_lm, exp_svm])
concat = concat.loc[~concat['model'].str.contains('PVC')]

for model in spearman_lm.columns:
    concat.loc[concat['model'].eq(model), 'Spearman correlation'] = spearman_lm[model].values

for model in spearman_svm.columns:
    concat.loc[concat['model'].eq(model), 'Spearman correlation'] = spearman_svm[model].values

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

# pairwise comparisons
pairs = [
    ('All binary', 'Binary A'),
    ('All binary', 'Binary T'),
    ('All binary', 'Binary N'),
    ('All categorical', 'Categorical A'),
    ('All categorical', 'Categorical T'),
    ('All categorical', 'Categorical N'),
    ('All continuous', 'Continuous A'),
    ('All continuous', 'Continuous T'),
    ('All continuous', 'Continuous N'),
    ('ATN SVM', 'Amyloid SVM'),
    ('ATN SVM', 'Tau SVM'),
    ('ATN SVM', 'GM SVM'),
    ('Amyloid SVM', 'Continuous A'),
    ('Tau SVM', 'Continuous T'),
    ('GM SVM', 'Continuous N'),
    ('ATN SVM', 'All continuous'),
    ]

positions = [1., 1.05, 1.1] * 4 + [1.2, 1.25, 1.3]
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='Baseline',
                             stats_vs_baseline=True, palette=palette,
                             n_train=n_train, n_test=n_test, font_file='arial.ttf',
                             stats_pairs=pairs, stats_pairs_positions=positions,
                             pivot_values = 'Spearman correlation', figsize=(7, 3.5))
set_labels_baseline_exp(fig)

plot_path = os.path.join('figures', 'spearman_correlations', 'results_combo_vs_baseline.svg')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)

#%% Plot spearman correlations (Combo vs. binary)

exp_lm = load_results('exp3_comboVSbinary_MMSE', 'results.csv')
exp_svm = load_results('exp1_svm_MMSE', 'results.csv')
spearman_lm = pd.read_csv(os.path.join('figures', 'spearman_correlations', 'combo_vs_binary.csv'))
spearman_svm = pd.read_csv(os.path.join('figures', 'spearman_correlations', 'svm.csv'))

concat = pd.concat([exp_lm, exp_svm])
concat = concat.loc[~concat['model'].str.contains('PVC')]

for model in spearman_lm.columns:
    concat.loc[concat['model'].eq(model), 'Spearman correlation'] = spearman_lm[model].values

for model in spearman_svm.columns:
    concat.loc[concat['model'].eq(model), 'Spearman correlation'] = spearman_svm[model].values

palette = (['#F7A934'] +
    ['#E899EE'] * 3 +
    ['#B74CBF'] +
    ['#F29D9D'] * 3 +
    ['#FC4646'] +
    ['#A5DBF2'] * 3 +
    ['#08A3E5'])

# plot
# plt.rcParams['font.size'] = 10
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='All binary',
                             stats_vs_baseline=True, palette=palette,
                             n_train=n_train, n_test=n_test, font_file='arial.ttf',
                             pivot_values='Spearman correlation', figsize=(7, 3.5))
set_labels_binary_exp(fig)

plot_path = os.path.join('figures', 'spearman_correlations', 'results_combo_vs_binary.svg')
plt.tight_layout()
fig.savefig(plot_path, dpi=300)
