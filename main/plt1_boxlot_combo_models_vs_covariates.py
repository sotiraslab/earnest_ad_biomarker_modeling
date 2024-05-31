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
exp2 = load_results('exp2_combo_atn_models_global_cognition', 'results.csv')
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')

# concatenate
concat = pd.concat([exp2, exp1])
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

positions = [0.75, 0.725, 0.7] * 4

# plot
# plt.rcParams['font.size'] = 10
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='Baseline',
                         stats_vs_baseline=True, palette=palette,
                         n_train=n_train, n_test=n_test, font_file='arial.ttf',
                         stats_pairs=pairs, stats_pairs_positions=positions)
set_labels_baseline_exp(fig)

plt.tight_layout()
fig.savefig(plot_path, dpi=300)
    