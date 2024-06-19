#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
from collections import Counter

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from atn_modeling import atn_predictor_classes
from atn_modeling.atn_predictor_instances import ATN_PREDICTORS
from common import load_results, set_font_properties

# ---- load selected models

models = load_results('exp2_combo_atn_models_global_cognition', 'models.pickle')

# ---- collect stuff for plotting

amyloid_cmap = LinearSegmentedColormap.from_list('Amyloid', colors=['#f7ebf1', '#882255'])

model_types = ['Binary A', 'Categorical A', 'Continuous A',
                'Binary T', 'Categorical T', 'Continuous T',
                'Binary N', 'Categorical N', 'Continuous N']
cmaps = [amyloid_cmap] * 3 + ['Greens'] * 3 + ['Blues'] * 3

# ---- function for converting "model_types" entry to something in ATN_PREDICTORS_DICT

def get_all_tested_model_names(model_type):
    measure, biomarker = model_type.split()
    measure = measure.lower()
    biomarker = {'A': 'amyloid', 'T': 'tau', 'N': 'neurodegeneration'}[biomarker]
    models = ATN_PREDICTORS[biomarker][measure]
    return [m.nickname for m in models]

# ---- collect data to plot
plot_dataframes = []

for model_type in model_types:
    data = []
    selected_models = [m[0].nickname for m in models[model_type]]
    counts = Counter(selected_models)
    data.append(counts)

    df = pd.DataFrame(data).T.sort_index()
    df = df.reindex(get_all_tested_model_names(model_type)).fillna(0).sort_values(0, ascending=False)
    df.index = df.index.str.replace('Amyloid ', '')
    plot_dataframes.append(df)

# ---- plot

# set font
set_font_properties('arial.ttf')

figure, axes = plt.subplots(nrows=3,
                            ncols=3,
                            figsize=(20, 14),
                            sharex=True,
                            sharey=True,
                            dpi=150)
# plt.subplots_adjust(wspace=2)

axes = axes.flatten()

for r, df in enumerate(plot_dataframes):
    ax = axes[r]
    cmap = plt.get_cmap(cmaps[r]).reversed()
    n_models_tested = len(df)
    n_models_used = df.iloc[:, 0].ne(0).sum()

    colors = np.zeros((n_models_tested, 4))
    colors[:, :] = [.64, .64, .64, 1.]
    colors[:n_models_used] = cmap((np.arange(n_models_used)/n_models_used) * .8)

    ax.pie(df[0], labels=None, colors=colors)


    # format
    ax.legend(bbox_to_anchor=(.9,.5), loc="center left", fontsize=14,
              labels = df.index)

bigfont=24
axes[0].set_ylabel('Amyloid', fontsize=bigfont)
axes[3].set_ylabel('Tau', fontsize=bigfont)
axes[6].set_ylabel('Neurodegeneration', fontsize=bigfont)

axes[0].set_title('Binary', fontsize=bigfont)
axes[1].set_title('Categorical', fontsize=bigfont)
axes[2].set_title('Continuous', fontsize=bigfont)

# save
plt.tight_layout()
plt.savefig('figures/feature_selection_pie_charts.svg', dpi=300)
