#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd

from common import load_results
from atn_modeling.helpers import results_boxplot

# load results
exp0 = load_results('exp0_individual_atn_models_global_cognition', 'results.csv')
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')
    
# concatenate data
exp1['variable_type'] = 'SVM'
exp1['name'] = exp1['model'].copy()
concat = pd.concat([exp0, exp1])
concat = concat[concat['biomarker'].isin([None, 'amyloid', 'tau', 'neurodegeneration']) | concat['name'].eq('Baseline')]

# plotting parameters
plot_path = os.path.join('figures', os.path.splitext(os.path.basename(__file__))[0] + '.svg')

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
                         palette=group['color'], order=order, font_file='arial.ttf')
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

# save
plt.tight_layout()
fig.savefig(plot_path, dpi=300)
