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
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')
exp10b = load_results('exp10b_test_all_predictors_with_csf_vs_baseline', 'results.csv')

# concatenate data
concat = pd.concat([exp10b, exp1])
concat = concat.loc[concat['model'].str.contains('Imaging') |
                    concat['model'].str.contains('CSF') |
                    concat['model'].eq('Baseline')]


# general resources
palette = (['gray'] +
    ['#F7A934'] * 2+
    ['#ba635d'] * 2+
    ['#7DA1D8'] * 2+
    ['#99C494'])
hatch = [False] + [False, True] * 3 + [False]
pairs = [('All binary [Imaging]', 'All binary [CSF]'),
         ('All categorical [Imaging]', 'All categorical [CSF]'),
         ('All continuous [Imaging]', 'All continuous [CSF]')]

# plot
n_train = concat['ntrain'].values[0]
n_test = concat['ntest'].values[0]
fig, _ = results_boxplot(concat, groupby='model', baseline='Baseline',
                          stats_vs_baseline=True, n_train=n_train, n_test=n_test,
                          font_file='arial.ttf', palette=palette, hatch=hatch,
                          stats_pairs=pairs)
set_labels_combo_vs_csf(fig)

plt.tight_layout()
fig.savefig('figures/combo_models_vs_basline_with_csf.svg', dpi=300)
