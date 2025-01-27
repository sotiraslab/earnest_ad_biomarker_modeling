#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import matplotlib.pyplot as plt
import pandas as pd

from atn_modeling import atn_predictor_classes
from atn_modeling.helpers import results_boxplot
from common import load_results, set_labels_pvc_combo

# load results
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')
exp2 = load_results('exp2_combo_atn_models_global_cognition', 'results.csv')

# concatenate data
concat = pd.concat([exp2, exp1])
keep = ['Baseline',
        'Binary T',
        'Binary T [PVC]',
        'All binary',
        'All binary [PVC]',
        'Categorical T',
        'Categorical T [PVC]',
        'All categorical',
        'All categorical [PVC]',
        'Continuous T',
        'Continuous T [PVC]',
        'All continuous',
        'All continuous [PVC]',
        'Tau SVM',
        'Tau SVM [PVC]',
        'ATN SVM',
        'ATN SVM [PVC]']
concat = concat.loc[concat['model'].isin(keep)]

# general resources
palette = (['gray'] +
    ['#FFDDAA'] * 2 +
    ['#F7A934'] * 2+
    ['#E899EE'] * 2 +
    ['#B74CBF'] * 2+
    ['#F29D9D'] * 2 +
    ['#FC4646'] * 2+
    ['#A5DBF2'] * 2 +
    ['#08A3E5'] * 2)
hatch = [False] + [False, True] * 8
pairs = list(zip(keep[1::2], keep[2::2]))

# plot
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, _ = results_boxplot(concat, groupby='model', baseline='Baseline',
                         stats_vs_baseline=True, palette=palette, hatch=hatch,
                         n_train=n_train, n_test=n_test, stats_pairs=pairs)
set_labels_pvc_combo(fig)

plt.tight_layout()
fig.savefig('figures/taupvc_combination_models.svg', dpi=300)
