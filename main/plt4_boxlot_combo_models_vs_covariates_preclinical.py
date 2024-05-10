#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import matplotlib.pyplot as plt
import pandas as pd

from common import load_results
from atn_modeling.helpers import results_boxplot

# load data
exp9a = load_results('exp9a_preclinical_combo_atn_models', 'results.csv')
exp9c = load_results('exp9c_preclinical_svms', 'results.csv')

# concatenate
concat = pd.concat([exp9a, exp9c])
concat = concat.loc[~concat['model'].str.contains('PVC')]

# output
plot_path = os.path.join('figures', os.path.splitext(os.path.basename(__file__))[0] + '.svg')

# general resources
palette = (['gray'] +
    ['#FFDDAA'] * 3 +
    ['#F7A934'] +
    ['#E59C9C'] * 3 +
    ['#ba635d'] +
    ['#AAC3E9'] * 3 +
    ['#7DA1D8'] +
    ['#B3E2AD'] * 3 +
    ['#99C494'])

# plot
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, _ = results_boxplot(concat, groupby='model', baseline='Baseline',
                         stats_vs_baseline=True, palette=palette,
                         n_train=n_train, n_test=n_test, font_file='arial.ttf')

plt.tight_layout()
fig.savefig(plot_path, dpi=300)
