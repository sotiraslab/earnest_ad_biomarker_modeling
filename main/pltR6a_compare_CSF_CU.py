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
from common import load_results, set_labels_combo_vs_csf

# load results
csf = load_results('expR6a_CSFcompare_MMSE_CU', 'results.csv')

# concatenate data
csf = csf.loc[csf['model'].str.contains('Imaging') |
              csf['model'].str.contains('CSF') |
              csf['model'].eq('Baseline')]


# general resources
palette = (['gray'] +
    ['#F7A934'] * 2+
    ['#B74CBF'] * 2+
    ['#FC4646'] * 2)
hatch = [False] + [False, True] * 3 + [False]
pairs = [('All binary [Imaging]', 'All binary [CSF]'),
         ('All categorical [Imaging]', 'All categorical [CSF]'),
         ('All continuous [Imaging]', 'All continuous [CSF]')]

# plot
n_train = csf['ntrain'].values[0]
n_test = csf['ntest'].values[0]
fig, stats = results_boxplot(csf, groupby='model', baseline='Baseline',
                             stats_vs_baseline=True, n_train=n_train, n_test=n_test,
                             font_file='arial.ttf', palette=palette, hatch=hatch,
                             stats_pairs=pairs, figsize=(3.5, 5))
set_labels_combo_vs_csf(fig)

# update text
fig.gca().set_ylabel('Prediction Error for ΔMMSE (RMSE)')

plt.tight_layout()
fig.savefig('figures/combo_models_vs_basline_with_csf_CU.svg', dpi=300)
