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
biomarkers = load_results('exp5b_comboVSsolo_MMSE_CU', 'results.csv')
svms = load_results('exp5a_svm_MMSE_CU', 'results.csv')

# concatenate
concat = pd.concat([biomarkers, svms])
concat = concat.loc[~concat['model'].str.contains('PVC')]

# output
plot_path = os.path.join('figures', os.path.splitext(os.path.basename(__file__))[0] + '.svg')

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
    ]

positions = [.6, .6, .64] * 4

# plot
# plt.rcParams['font.size'] = 10
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='Baseline',
                         stats_vs_baseline=True, palette=palette,
                         n_train=n_train, n_test=n_test, font_file='arial.ttf',
                         stats_pairs=pairs, stats_pairs_positions=positions,
                         figsize=(7, 3.5))
set_labels_baseline_exp(fig)

# update text
fig.gca().set_ylabel('Prediction Error for ΔMMSE (RMSE)')

plt.tight_layout()
fig.savefig(plot_path, dpi=300)
