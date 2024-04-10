#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os
import pickle

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from pingouin import ttest
from statsmodels.stats.multitest import multipletests

import sys
sys.path.append('../1_modeling')

# ---- locate output

if not os.path.exists('../../outputs/additional_plots'):
    os.mkdir('../../outputs/additional_plots')

# ---- load saved models from binary combination experiment

separate_models_path = '../../outputs/exp0_individual_atn_models_global_cognition/models.pickle'
combo_models_path = '../../outputs/exp2_combo_atn_models_global_cognition/models.pickle'
with open(separate_models_path, 'rb') as f:
    sep_models = pickle.load(f)
    
with open(combo_models_path, 'rb') as f:
    com_models = pickle.load(f)
    
sep_results = pd.read_csv('../../outputs/exp0_individual_atn_models_global_cognition/results.csv')
com_results = pd.read_csv('../../outputs/exp2_combo_atn_models_global_cognition/results.csv')

# ---- select all tau models

pvc_model_names = []
nopvc_model_names = []

for k, v in sep_models.items():
    m = v[0]
    if m.atn == 'tau':
        nopvc_model_names.append(k)
    elif m.atn == 'taupvc':
        pvc_model_names.append(k)
    else:
        pass
    
# ---- plot: for separate models

from helpers import results_boxplot_pairwise

# set styling
colors = {None: 'Gray',
          'tau': '#117733',
          'taupvc': '#b48dc4'}
vartypes = {'binary': "BIN",
            'categorical': "CAT",
            'continuous': 'CONT'}

# order PVC and non PVC next to each other
interleave = [val for pair in zip(nopvc_model_names, pvc_model_names) for val in pair]
interleave = ['Baseline'] + interleave

# comparisons for PVC and non-PVC
pairs = list(zip(nopvc_model_names, pvc_model_names))

# construct plot data
plot_data = sep_results.copy()
plot_data = plot_data.loc[plot_data['name'].isin(nopvc_model_names) | 
                          plot_data['name'].isin(pvc_model_names) |
                          plot_data['name'].eq('Baseline')]
group = plot_data.groupby('name').agg({'name': 'first', 'biomarker': 'first',
                                       'variable_type': 'first', 'rmse': 'mean'})
group.index.name = 'Index'
group['color'] = group['biomarker'].map(colors)
group = group.reindex(interleave)

# run!
fig, sep_stats = results_boxplot_pairwise(plot_data, groupby='name', baseline='Baseline',
                                          palette=group['color'],
                                          order=group.index,
                                          n_train=plot_data['ntrain'].values[-1],
                                          n_test=plot_data['ntest'].values[-1],
                                          pairs=pairs)

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

# ---- plot: for combination models