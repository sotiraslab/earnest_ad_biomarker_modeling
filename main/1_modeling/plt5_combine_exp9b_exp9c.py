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

def main():
    
    this_dir = os.path.dirname(os.path.abspath(__file__))
    expB_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp9b_preclinical_combo_atn_models_vs_binary'))
    expC_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp9c_preclinical_svms'))
    
    expB_results = os.path.join(expB_folder, 'results.csv')
    expC_results = os.path.join(expC_folder, 'results.csv')
    
    if not os.path.exists(expB_results):
        raise RuntimeError('Results for exp9b are missing.')
        
    if not os.path.exists(expC_results):
        short_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp9c_preclinical_svms_short'))
        short_results = os.path.join(short_folder, 'results.csv')
        
        if os.path.exists(short_results):
            warnings.warn('The full results are missing for exp9c (SVMs), but the short results are present.  Using the latter!',
                          RuntimeWarning)
            expC_results = short_results
        else:
            raise RuntimeError('Results for exp9c are missing.')
        
    # concatenate data
    expB = pd.read_csv(expB_results)
    expC = pd.read_csv(expC_results)
    concat = pd.concat([expB, expC])
    concat = concat.loc[~concat['model'].str.contains('PVC')]
    
    # output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    # general resources
    plot_path = os.path.join(output_folder, 'boxplot.svg')
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
    
    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
    
if __name__ == '__main__':
    main()