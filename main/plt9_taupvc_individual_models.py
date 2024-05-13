#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from atn_modeling import atn_predictor_classes
from atn_modeling.helpers import results_boxplot
from common import load_results

# ------ load results ------

exp0 = load_results('exp0_individual_atn_models_global_cognition', 'results.csv')
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')

# ---- select all tau models

pvc_model_names = list(exp0.loc[exp0['biomarker'] == 'taupvc', 'name'].unique())
nopvc_model_names = list(exp0.loc[exp0['biomarker'] == 'tau', 'name'].unique())

pvc_model_names += ['Tau SVM [PVC]']
nopvc_model_names += ['Tau SVM']

# ---- plot: for separate models

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

# APPLYING A RANGE LIMIT ON THE GRAPH

# For this seed, there is one iteration (repeat 8, fold 7)
# where the accuracy for MTT (UZ) [PVC] is >0.89,
# much worse than the rest.  It messes up the range of the graph.

# This bit adds a manual xlim, and filters that point.
# It should report if a point is actually removed.

# Note that this doesn't affect the box plot computation,
# nor the statistics.

# Make a note of this in the figure caption if included in published content.
limit = plot_data['rmse'].mean() + (plot_data['rmse'].std() * 4)
omit_mask = plot_data['rmse'] >= limit
omitted = omit_mask.sum()
omitted_models = list(plot_data.loc[omit_mask, 'name'].unique())
if omitted:
    print('!!! !!! !!!')
    print('NOTE: points omitted from figure for visualization purposes.')
    print(f'  - number omitted: {omitted}')
    print(f'  - models affected: {omitted_models}')
    print('!!! !!! !!!')

ax.set_xlim(None, limit)

# save
plt.tight_layout()
fig.savefig('figures/taupvc_individual_models.svg')
