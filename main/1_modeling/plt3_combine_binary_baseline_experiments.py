#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import warnings

import matplotlib.pyplot as plt
import pandas as pd

from helpers import results_boxplot

# constats
this_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = os.path.join(this_dir, '..', '..', 'outputs')

# more helper functions
def contains_results(experiment):
    return os.path.exists(os.path.join(output_dir, experiment, 'results.csv'))

def locate_outfolder(experiment):
    outputs = os.listdir(output_dir)
    long_outputs = [s for s in outputs if '_short' not in s]
    search1 = [s for s in long_outputs if experiment in s and contains_results(s)]
    search2 = [s for s in outputs if experiment in s and contains_results(s)]
    if search1:
        return os.path.join(output_dir, search1[0])
    elif search2:
        warnings.warn(f'Only found short results for "{experiment}".  Using those...')
        return os.path.join(output_dir, search2[0])
    else:
        raise ValueError(f'Cannot find any results for "{experiment}".')
        
experiment_keys = [
    'memory',
    'executive_functioning',
    'language',
    'visuospatial',
    'longitudinal_cognition'
    ]

for key in experiment_keys:
    try:
        lm_directory = locate_outfolder(f'{key}_vs_binary')
        svm_directory = locate_outfolder(f'svms_{key}')
    except ValueError:
        continue
    
    # concatenate data
    results_lm = pd.read_csv(os.path.join(lm_directory, 'results.csv'))
    results_svm = pd.read_csv(os.path.join(svm_directory, 'results.csv'))
    concat = pd.concat([results_lm, results_svm])
    concat = concat.loc[~concat['model'].str.contains('PVC')]
    
    # output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)

    # general resources
    plot_path = os.path.join(output_folder, f'{key}.svg')
    palette = (['gray'] +
        ['#E59C9C'] * 3 +
        ['#ba635d'] +
        ['#AAC3E9'] * 3 +
        ['#7DA1D8'] +
        ['#B3E2AD'] * 3 +
        ['#99C494'])

    # plot
    n_train = concat['ntrain'].values[0]
    n_test = concat['ntest'].values[0]
    fig, _ = results_boxplot(concat, groupby='model', baseline='All binary',
                              stats_vs_baseline=True, palette=palette,
                              n_train=n_train, n_test=n_test)
    plt.title(key)

    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)


