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
exp12 = load_results('exp10a_test_all_predictors_with_csf', 'results.csv')
exp1 = load_results('exp1_svms_global_cognition', 'results.csv')

# concatenate data
exp1['variable_type'] = 'SVM'
exp1['name'] = exp1['model'].copy()
concat = pd.concat([exp12, exp1])
bmarkers = [None, 'amyloid', 'tau', 'neurodegeneration', 'csf_tau', 'csf_amyloid', 'csf_neurodegeneration']
concat = concat[concat['biomarker'].isin(bmarkers) | concat['name'].eq('Baseline')]

colors = {'amyloid': '#882255',
          'tau': '#117733',
          'neurodegeneration': '#332288',
          'csf_amyloid': '#882255',
          'csf_tau': '#117733',
          'csf_neurodegeneration': '#332288',
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
group['csf'] = group['biomarker'].str.contains('csf')
group['csf'] = group['csf'].fillna(False)
group = group.loc[order, :]
group['sorter'] = group['rmse']
group.loc[group['name'].eq('Baseline'), 'sorter'] += 100
group = group.sort_values('sorter', ascending=False)
order = group['name']

fig, _ = results_boxplot(data, groupby='name', baseline='Baseline',
                         palette=group['color'], order=order, font_file='arial.ttf',
                         hatch=group['csf'])
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
    mpatches.Patch(color=colors['amyloid'], label='Amyloid (Imaging)'),
    mpatches.Patch(facecolor=colors['amyloid'], edgecolor='white', label='Amyloid (CSF)', hatch='//'),
    mpatches.Patch(color=colors['tau'], label='Tau (Imaging)'),
    mpatches.Patch(facecolor=colors['tau'], edgecolor='white', label='Tau (CSF)', hatch='//'),
    mpatches.Patch(color=colors['neurodegeneration'], label='ND (Imaging)'),
    mpatches.Patch(facecolor=colors['neurodegeneration'], edgecolor='white', label='ND (CSF)', hatch='//'),
    ],
    loc='lower right',
    bbox_to_anchor=(1, 1),
    ncol=3,
    frameon=False)

# save
plt.tight_layout()
fig.savefig('figures/individual_models_vs_basline_with_csf.svg', dpi=300)
