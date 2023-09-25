#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
from collections import Counter
import os
import pickle 
import sys

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# must reload the model classes in order for pickle to read them
sys.path.append('../1_modeling')

# ---- load selected models
def load_results(folder):
    ppath = os.path.join(folder, 'models.pickle')
    with open(ppath, 'rb') as f:
        models = pickle.load(f)
        
    return models

models_all = load_results('../1_modeling/control_baseline')
models_cn = load_results('../1_modeling/control_baseline_CN')
models_dem = load_results('../1_modeling/control_baseline_DEM')
models_long = load_results('../1_modeling/control_baseline_longitudinal')

# ---- collect stuff for plotting
experiments = [models_all, models_cn, models_dem, models_long]
experiment_names = ['All', 'CDR=0', 'CDR>0', 'All (longitudinal)']

model_types = ['Binary A', 'Categorical A', 'Continuous A',
               'Binary T', 'Categorical T', 'Continuous T',
               'Binary N', 'Categorical N', 'Continuous N']
cmaps = ['Purples'] * 3 + ['Blues'] * 3 + ['Greens'] * 3


# ---- collect data to plot
plot_dataframes = []

for model_type in model_types:
    data = []
    for experiment in experiments:
        selected_models = [str(m[0]) for m in experiment[model_type]]
        counts = Counter(selected_models)
        data.append(counts)
        
    df = pd.DataFrame(data, index=experiment_names).T.sort_index()
    df = df.fillna(0)
    df/= df.sum()
    plot_dataframes.append(df)

# ---- plot

# font
font_prop = fm.FontProperties(fname='../../fonts/arial.ttf')
plt.rcParams.update({
    'font.family': font_prop.get_name(),
    'font.size': 14})

figure, axes = plt.subplots(nrows=3,
                            ncols=3,
                            figsize=(20, 12),
                            sharex=True,
                            sharey=True,
                            dpi=300)
axes = axes.flatten()

cols = len(experiment_names)
for r, df in enumerate(plot_dataframes):
    ax = axes[r]
    unique_bars = len(df) 
    bottom = np.zeros(cols)
    cmap = plt.get_cmap(cmaps[r]).reversed()
    colors = cmap((np.arange(unique_bars)/unique_bars) * .8) 
    for i in range(unique_bars):
        
        # TODO: get better labels/short names for models
        label = df.index[i][:20]
        
        height = df.iloc[i, :]
        ax.bar(x=range(cols), height=height, bottom=bottom, label=label, color=colors[i, :])
        bottom += height
        
    # format
    ax.set_ylim(0,1)
    ax.set_yticks([0, .5, 1])
    ax.legend(bbox_to_anchor=(1,1), loc="upper left", fontsize=12)
    ax.set_xticks(range(cols))
    ax.set_xticklabels(experiment_names, rotation=45, ha='right')

bigfont=24
axes[0].set_ylabel('Amyloid', fontsize=bigfont)
axes[3].set_ylabel('Tau', fontsize=bigfont)
axes[6].set_ylabel('Neurodegeneration', fontsize=bigfont)

axes[0].set_title('Binary', fontsize=bigfont)
axes[1].set_title('Categorical', fontsize=bigfont)
axes[2].set_title('Continuous', fontsize=bigfont)
    
plt.tight_layout()
plt.savefig('selected_models.png', dpi=300)