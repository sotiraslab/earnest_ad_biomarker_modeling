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
csf = load_results('exp7a_CSFshow_MMSE', 'results.csv')
svm = load_results('exp1_svm_MMSE', 'results.csv')

# concatenate data
svm['variable_type'] = 'SVM'
svm['name'] = svm['model'].copy()
concat = pd.concat([csf, svm])
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
varcolors = {'binary': "#F7A934",
            'categorical': "#B74CBF",
            'continuous': '#FC4646',
            'SVM': '#08A3E5'}

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

# make text a little smaller for this plot
plt.rcParams.update({'font.size': 9})

fig, _ = results_boxplot(data, groupby='name', baseline='Baseline',
                         palette=group['color'], order=order, font_file='arial.ttf',
                         hatch=group['csf'])
ax = fig.axes[0]
for i, var in enumerate(group['variable_type']):
    if var is None:
        continue
    label = vartypes[var]
    color = varcolors[var]
    ax.text(i, 1.02, label, ha='center', va='bottom', color=color, rotation=90,
            transform=ax.get_xaxis_transform())

# manually make legend
ax.legend(handles = [
    mpatches.Patch(color=colors['amyloid'], label='Amyloid (Imaging)'),
    mpatches.Patch(facecolor=colors['amyloid'], edgecolor='white', label='Amyloid (CSF)', hatch='//'),
    mpatches.Patch(color=colors['tau'], label='Tau (Imaging)'),
    mpatches.Patch(facecolor=colors['tau'], edgecolor='white', label='Tau (CSF)', hatch='//'),
    mpatches.Patch(color=colors['neurodegeneration'], label='ND (Imaging)'),
    mpatches.Patch(facecolor=colors['neurodegeneration'], edgecolor='white', label='ND (CSF)', hatch='//'),
    ],
    loc='lower center',
    bbox_to_anchor=(0.5, 1.07),
    ncol=3,
    frameon=False)

# update text
ax.set_ylabel('Prediction Error for ΔMMSE (RMSE)')

# save
plt.tight_layout()
fig.savefig('figures/individual_models_vs_basline_with_csf.svg', dpi=300)
