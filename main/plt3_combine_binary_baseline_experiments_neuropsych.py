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
from common import locate_outfolder, set_labels_binary_exp

output_dir = 'outputs'

# more helper functions


experiment_keys = [
    'memory',
    'executive_functioning',
    'language',
    'visuospatial'
    ]


stats_dict = {}
for key in experiment_keys:
    try:
        lm_directory = locate_outfolder(f'{key}_vs_binary')
        svm_directory = locate_outfolder(f'svms_{key}')
    except ValueError:
        print(f'NO RESULTS FOUND FOR KEY: {key}')
        continue

    # concatenate data
    results_lm = pd.read_csv(os.path.join(lm_directory, 'results.csv'))
    results_svm = pd.read_csv(os.path.join(svm_directory, 'results.csv'))
    concat = pd.concat([results_lm, results_svm])
    concat = concat.loc[~concat['model'].str.contains('PVC')]

    # output
    plot_path = os.path.join('figures', f'plt3_boxplot_vs_binary_{key}.svg')

    # general resources
    palette = (
        ['#F7A934'] +
        ['#E899EE'] * 3 +
        ['#B74CBF'] +
        ['#F29D9D'] * 3 +
        ['#FC4646'] +
        ['#A5DBF2'] * 3 +
        ['#08A3E5']
        )

    # plot
    n_train = concat['ntrain'].values[0]
    n_test = concat['ntest'].values[0]
    fig, stats = results_boxplot(concat, groupby='model', baseline='All binary',
                                 stats_vs_baseline=True, palette=palette,
                                 n_train=n_train, n_test=n_test)
    set_labels_binary_exp(fig)
    stats_dict[key] = stats['baseline']
    plt.title(key)

    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)


