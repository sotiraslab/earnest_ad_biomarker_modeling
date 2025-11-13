#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import matplotlib.pyplot as plt
import pandas as pd

from common import load_results, locate_outfolder, set_labels_baseline_exp
from atn_modeling.helpers import results_boxplot

experiment_keys = [
    'MEM',
    'EXF',
    'LAN',
    'VSP',
    'PC1'
    ]

stats_dict = {}
for key in experiment_keys:
    try:
        lm_directory = locate_outfolder(f'comboVSsolo_DeltaPHC_{key}')
        svm_directory = locate_outfolder(f'svm_{key}')
    except ValueError:
        print(f'NO RESULTS FOUND FOR KEY: {key}')
        continue

    # load data
    results_lm = pd.read_csv(os.path.join(lm_directory, 'results.csv'))
    results_svm = pd.read_csv(os.path.join(svm_directory, 'results.csv'))
    concat = pd.concat([results_lm, results_svm])
    concat = concat.loc[~concat['model'].str.contains('PVC')]

    # output
    plot_path = os.path.join('figures', f'plt3_boxplot_vs_solo_{key}.svg')


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

    # plot
    # plt.rcParams['font.size'] = 10
    n_train = concat['ntrain'].values[0]
    n_test = concat['ntest'].values[0]
    fig, stats = results_boxplot(concat, groupby='model', baseline='Baseline',
                             stats_vs_baseline=True, palette=palette,
                             n_train=n_train, n_test=n_test, font_file='arial.ttf',
                             stats_pairs=pairs)
    set_labels_baseline_exp(fig)

    # update text
    fig.gca().set_ylabel(f'Prediction Error for ΔPHC-{key} (RMSE)')

    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
