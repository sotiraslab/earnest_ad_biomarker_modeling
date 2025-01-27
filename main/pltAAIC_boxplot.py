#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 09:54:56 2024

@author: tom.earnest
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_results
from atn_modeling.helpers import results_boxplot

# load data
data = load_results('expAAIC_svms', 'results.csv')
    
# output
plot_path = os.path.join('figures', os.path.splitext(os.path.basename(__file__))[0] + '.svg')

# colors
palette = [
    'Gray',
    
    '#fe6100',
    '#dc267f',
    '#648fff',
    '#785ef0',
    
    '#fe6100',
    '#dc267f',
    '#648fff',
    '#785ef0',]

# plot
# plt.rcParams['font.size'] = 10
n_train = data['ntrain'].values[0]
n_test = data['ntest'].values[0]
fig, stats = results_boxplot(data, groupby='model', baseline='Baseline',
                             stats_vs_baseline=False, palette=palette,
                             n_train=n_train, n_test=n_test, font_file='arial.ttf')
ax = plt.gca()
ax.set_yticklabels(
    ['Age+Sex+APOE4',
     '+ Amyloid SUVR',
     '+ Tau SUVR',
     '+ Temporal volume',
     '+ ATN composites',
     '+ Amyloid ROIs',
     '+ Tau ROIs',
     '+ GM ROIs',
     '+ ATN ROIs'])

ax.axhspan(-0.5, 3.5, alpha=0.3, color='#ffd966')
ax.axhspan(3.5, 7.5, alpha=0.3, color='#6fa8dc')
ax.text(0.7, 3.75,
        'Composites',
        ha='right',
        va='center',
        color='#0f7bdd',
        fontsize=12)

ax.text(0.7, -0.25,
        'Regional',
        ha='right',
        va='center',
        color='#ce7e00',
        fontsize=12)

plt.tight_layout()
fig.savefig(plot_path, dpi=300)

# recreate plot to get stats
pairs = [
    # biomarkers vs each other
    ('Amyloid biomarker', 'Tau biomarker'),
    ('Amyloid biomarker', 'GM biomarker'),
    ('Tau biomarker', 'GM biomarker'),
    
    # unimodal biomarkers vs multimodal
    ('ATN biomarker', 'Amyloid biomarker'),
    ('ATN biomarker', 'Tau biomarker'),
    ('ATN biomarker', 'GM biomarker'),
    
    # ROIs vs each other
    ('Amyloid ROI', 'Tau ROI'),
    ('Amyloid ROI', 'GM ROI'),
    ('Tau ROI', 'GM ROI'),
    
    # unimodal ROIs vs multimodal ROIs
    ('ATN ROI', 'Amyloid ROI'),
    ('ATN ROI', 'Tau ROI'),
    ('ATN ROI', 'GM ROI'),
    
    # biomarkers vs ROIs
    ('Amyloid biomarker', 'Amyloid ROI'),
    ('Tau biomarker', 'Tau ROI'),
    ('GM biomarker', 'GM ROI'),
    ('ATN biomarker', 'ATN ROI')
    ]


_, stats = results_boxplot(data, groupby='model', baseline='Baseline',
                           stats_vs_baseline=True, palette=palette,
                           n_train=n_train, n_test=n_test, font_file='arial.ttf',
                           stats_pairs=pairs)
stats_bl = stats['baseline']
stats_pairs = stats['pairs']

# nice version of table
nice = stats_pairs[['a', 'b', 'mean_a', 'mean_b', 't', 'p-val-fdr']].copy()
nice['mean_a'] = nice['mean_a'].round(2)
nice['mean_b'] = nice['mean_b'].round(2)
nice['t'] = nice['t'].round(2)
nice['p-val-fdr'] = nice['p-val-fdr'].round(2)
stars =  pd.cut(nice['p-val-fdr'], [-np.inf, 0.001, 0.01, 0.05, np.inf], labels=['***', '**', '*', ''])
nice['p-val-fdr'] = nice['p-val-fdr'].astype(str)
nice.loc[nice['p-val-fdr'].eq('0.0'), 'p-val-fdr'] = '<0.001'
nice['p-val-fdr'] = nice['p-val-fdr'] + ' ' + stars.astype(str)
nice['p-val-fdr'] = nice['p-val-fdr'].str.strip()

remap = {
    'Amyloid biomarker':'Amyloid SUVR',
    'Tau biomarker': 'Tau SUVR',
    'GM biomarker': 'Temporal volume',
    'ATN biomarker': 'ATN biomarkers',
    'Amyloid ROI': 'Amyloid ROIs',
    'Tau ROI': 'Tau ROIs',
    'GM ROI': 'GM ROIs',
    'ATN ROI': 'ATN ROIs'
    }
nice['a'] = nice['a'].map(remap)
nice['b'] = nice['b'].map(remap)
nice = nice[['a', 'mean_a', 'b', 'mean_b', 't', 'p-val-fdr']]