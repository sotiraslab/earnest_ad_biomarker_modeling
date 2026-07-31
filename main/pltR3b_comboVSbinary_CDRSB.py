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
from common import load_results, set_labels_binary_exp

# concatenate data
biomarkers = load_results('expR3b_comboVSbinary_CDRSB', 'results.csv')
svm = load_results('expR3c_svm_CDRSB', 'results.csv')
concat = pd.concat([biomarkers, svm])
concat = concat.loc[~concat['model'].str.contains('PVC')]

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

# plot
# plt.rcParams['font.size'] = 10
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='All binary',
                             stats_vs_baseline=True, palette=palette,
                             n_train=n_train, n_test=n_test, font_file='arial.ttf',
                             pivot_values='rmse', figsize=(7, 3.5))
set_labels_binary_exp(fig)
fig.gca().set_ylabel('Prediction Error for ΔCDR-SB (RMSE)')

plt.tight_layout()
fig.savefig(plot_path, dpi=300)
