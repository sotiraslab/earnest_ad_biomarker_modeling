#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os
import warnings

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

import sys
sys.path.append('../1_modeling')

# ---- locate output

if not os.path.exists('../../outputs/additional_plots'):
    os.mkdir('../../outputs/additional_plots')

# ------ load results ------

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
        
exp0 = pd.read_csv(exp0_results)
exp1 = pd.read_csv(exp1_results)

# ---- select all tau models

pvc_model_names = list(exp0.loc[exp0['biomarker'] == 'taupvc', 'name'].unique())
nopvc_model_names = list(exp0.loc[exp0['biomarker'] == 'tau', 'name'].unique())

pvc_model_names += ['Tau SVM [PVC]']
nopvc_model_names += ['Tau SVM']

# ---- plot: for separate models

from helpers import results_boxplot

# set styling
colors = {None: 'Gray',
          'tau': '#117733',
          'taupvc': '#b48dc4'}
vartypes = {'binary': "BIN",
            'categorical': "CAT",
            'continuous': 'CONT',
            'SVM': 'SVM'}

# order PVC and non PVC next to each other
interleave = [val for pair in zip(nopvc_model_names, pvc_model_names) for val in pair]
interleave = ['Baseline'] + interleave

# comparisons for PVC and non-PVC
pairs = list(zip(nopvc_model_names, pvc_model_names))

# concatenate data
exp0 = pd.read_csv(exp0_results)
exp1 = pd.read_csv(exp1_results)
exp1['variable_type'] = 'SVM'
exp1['name'] = exp1['model'].copy()
concat = pd.concat([exp0, exp1])
concat = concat[concat['name'].isin(interleave)]

plot_data = concat.copy()
plot_data = plot_data.loc[plot_data['name'].isin(nopvc_model_names) | 
                          plot_data['name'].isin(pvc_model_names) |
                          plot_data['name'].eq('Baseline')]
group = plot_data.groupby('name').agg({'name': 'first', 'biomarker': 'first',
                                       'variable_type': 'first', 'rmse': 'mean'})
group.index.name = 'Index'
group['color'] = group['biomarker'].map(colors)
group = group.reindex(interleave)

# run!
fig, stats = results_boxplot(plot_data, groupby='name', baseline='Baseline',
                             palette=group['color'],
                             order=group.index,
                             n_train=plot_data['ntrain'].values[-1],
                             n_test=plot_data['ntest'].values[-1],
                             stats_pairs=pairs)

# add vartypes to side
ax = fig.axes[0]
for i, var in enumerate(group['variable_type']):
    y = (len(group) - 1) - i
    if var is None:
        continue
    label = vartypes[var]
    ax.text(1.05, y, label, ha='left', va='center', transform=ax.get_yaxis_transform())

# make a legend
ax.legend(handles = [
    mpatches.Patch(color=colors['tau'], label='Tau'),
    mpatches.Patch(color=colors['taupvc'], label='Tau [PVC]'),
    ],
    loc='lower left',
    bbox_to_anchor=(0, 1),
    ncol=2,
    frameon=False)

# save
plt.tight_layout()
fig.savefig('../../outputs/additional_plots/tau_pvc_indvl_models.svg')
