#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd

from experiments import experiment_test_all_atn_predictors
from helpers import results_boxplot

def parse():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--noplot', action='store_true')
    
    return parser.parse_args()

def main(rerun=False, replot=True):
    
    # setup output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    this_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    # check existing output
    results_path = os.path.join(output_folder, 'results.csv')
    modelsave_path = os.path.join(output_folder, 'models.pickle')
    results_exist = os.path.exists(results_path)
        
    # run experiment
    if rerun or (not results_exist):
        dataset = pd.read_csv('../../outputs/maindata/maindata.csv')
        _ = experiment_test_all_atn_predictors(dataset, target='PHC_GLOBAL',
                                               covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                               stratify='CDRBinned',
                                               repeats=10, splits=10,
                                               seed=0,
                                               savepath=results_path,
                                               savemodels=modelsave_path)
        
    results = pd.read_csv(results_path)
        
    if replot:
        
        # general resources
        plot_path = os.path.join(output_folder, 'boxplot.svg')
        colors = {'amyloid': '#882255',
                  'tau': '#117733',
                  'neurodegeneration': '#332288',
                  None: 'Gray'}
        vartypes = {'binary': "BIN",
                    'categorical': "CAT",
                    'continuous': 'CONT'}
        
        # 1. grouped by biomarker
        data = results.loc[~results['name'].str.contains('PVC'), :]
        order = data['name'].unique()
        group = data.groupby('name').agg({'name': 'first', 'biomarker': 'first',
                                          'variable_type': 'first', 'rmse': 'mean'})
        group['color'] = group['biomarker'].map(colors)
        group = group.loc[order, :]
        
        fig, _ = results_boxplot(data, groupby='name', baseline='Baseline',
                                 stats=False, palette=group['color'])
        ax = fig.axes[0]
        for i, var in enumerate(group['variable_type']):
            y = (len(group) - 1) - i
            if var is None:
                continue
            label = vartypes[var]
            ax.text(1.05, y, label, ha='left', va='center', transform=ax.get_yaxis_transform())
        plt.tight_layout()
        fig.savefig(plot_path, dpi=300)

        # 2. ordered by accuracy
        plot_path = os.path.join(output_folder, 'boxplot_sorted.svg')
        group = group.sort_values('rmse', ascending=False)
        order = group['name']
        
        fig, _ = results_boxplot(data, groupby='name', baseline='Baseline',
                                 stats=False, palette=group['color'], order=order)
        ax = fig.axes[0]
        for i, var in enumerate(group['variable_type']):
            y = (len(group) - 1) - i
            if var is None:
                continue
            label = vartypes[var]
            ax.text(1.05, y, label, ha='left', va='center', transform=ax.get_yaxis_transform())
        plt.tight_layout()
        fig.savefig(plot_path, dpi=300)
    
if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    kwargs['replot'] = not kwargs['noplot']
    del kwargs['noplot']
    main(**kwargs)