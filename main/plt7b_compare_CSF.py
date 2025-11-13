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
svm = load_results('exp1_svm_MMSE', 'results.csv')
csf = load_results('exp7b_CSFcompare_MMSE', 'results.csv')

# concatenate data
concat = pd.concat([csf, svm])
concat = concat.loc[concat['model'].str.contains('Imaging') |
                    concat['model'].str.contains('CSF') |
                    concat['model'].eq('Baseline')]


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
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, stats = results_boxplot(concat, groupby='model', baseline='Baseline',
                             stats_vs_baseline=True, n_train=n_train, n_test=n_test,
                             font_file='arial.ttf', palette=palette, hatch=hatch,
                             stats_pairs=pairs)
set_labels_combo_vs_csf(fig)

plt.tight_layout()
fig.savefig('figures/combo_models_vs_basline_with_csf.svg', dpi=300)
