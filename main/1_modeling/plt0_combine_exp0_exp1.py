#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import warnings

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from helpers import results_boxplot

def main():
    
    this_dir = os.path.dirname(os.path.abspath(__file__))
    exp0_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp0_individual_atn_models_global_cognition'))
    exp1_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp1_svms_global_cognition'))
    
    exp0_results = os.path.join(exp0_folder, 'results.csv')
    exp1_results = os.path.join(exp1_folder, 'results.csv')
    
    if not os.path.exists(exp0_results):
        raise RuntimeError('Results for exp0 are missing.')
        
    if not os.path.exists(exp1_results):
        short_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp1_svms_global_cognition_short'))
        short_results = exp1_results = os.path.join(short_folder, 'results.csv')
        
        if os.path.exists(short_results):
            warnings.warn('The full results are missing for exp1 (SVMs), but the short results are present.  Using the latter!',
                          RuntimeWarning)
            exp1_results = short_results
        else:
            raise RuntimeError('Results for exp1 are missing.')
        
    # concatenate data
    exp0 = pd.read_csv(exp0_results)
    exp1 = pd.read_csv(exp1_results)
    exp1['variable_type'] = 'SVM'
    exp1['name'] = exp1['model'].copy()
    concat = pd.concat([exp0, exp1])
    concat = concat[concat['biomarker'].isin([None, 'amyloid', 'tau', 'neurodegeneration']) | concat['name'].eq('Baseline')]
        
    # output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
    
    # general resources
    plot_path = os.path.join(output_folder, 'boxplot.svg')
    colors = {'amyloid': '#882255',
              'tau': '#117733',
              'neurodegeneration': '#332288',
              None: 'Gray'}
    vartypes = {'binary': "BIN",
                'categorical': "CAT",
                'continuous': 'CONT',
                'SVM': 'SVM'}
    varcolors = {'binary': "#011f4b",
                'categorical': "#03396c",
                'continuous': '#005b96',
                'SVM': '#6497b1'}
    
    # group by model and plot
    data = concat.copy()
    order = data['name'].unique()
    group = data.groupby('name').agg({'name': 'first', 'biomarker': 'first',
                                      'variable_type': 'first', 'rmse': 'mean'})
    group['color'] = group['biomarker'].map(colors)
    group = group.loc[order, :]
    group = group.sort_values('rmse', ascending=False)
    order = group['name']
    
    fig, _ = results_boxplot(data, groupby='name', baseline='Baseline',
                             palette=group['color'], order=order)
    ax = fig.axes[0]
    for i, var in enumerate(group['variable_type']):
        y = (len(group) - 1) - i
        if var is None:
            continue
        label = vartypes[var]
        color = varcolors[var]
        ax.text(1.05, y, label, ha='left', va='center',
                transform=ax.get_yaxis_transform(), color=color)
        
    # manually make legend
    ax.legend(handles = [
        mpatches.Patch(color=colors['amyloid'], label='Amyloid'),
        mpatches.Patch(color=colors['tau'], label='Tau'),
        mpatches.Patch(color=colors['neurodegeneration'], label='Neurodegeneration'),
        ],
        loc='lower right',
        bbox_to_anchor=(1, 1),
        ncol=3,
        frameon=False)
    
    plt.tight_layout()
    fig.savefig(plot_path, dpi=300)
    
if __name__ == '__main__':
    main()