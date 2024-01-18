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

# also import the full list of models
from model_classes import ATN_PREDICTORS_DICT

# ---- load selected models
def load_results(folder):
    ppath = os.path.join(folder, 'models.pickle')
    with open(ppath, 'rb') as f:
        models = pickle.load(f)

    return models

models = load_results('../1_modeling/control_baseline')

# ---- collect stuff for plotting

model_types = ['Binary A', 'Categorical A', 'Continuous A',
                'Binary T', 'Categorical T', 'Continuous T',
                'Binary N', 'Categorical N', 'Continuous N']
cmaps = ['Purples'] * 3 + ['Blues'] * 3 + ['Greens'] * 3

# ---- function for converting "model_types" entry to something in ATN_PREDICTORS_DICT

def get_all_test_model_names(model_type):
    measure, biomarker = model_type.split()
    measure = measure.lower()
    biomarker = {'A': 'amyloid', 'T': 'tau', 'N': 'neurodegeneration'}[biomarker]
    models = ATN_PREDICTORS_DICT[biomarker][measure].values()
    return [m.nickname for m in models]


# ---- collect data to plot
plot_dataframes = []

for model_type in model_types:
    data = []
    selected_models = [m[0].nickname for m in models[model_type]]
    counts = Counter(selected_models)
    data.append(counts)

    df = pd.DataFrame(data).T.sort_index()
    df = df.reindex(get_all_test_model_names(model_type)).fillna(0).sort_values(0, ascending=False)
    df.index = df.index.str.replace('Amyloid ', '')
    plot_dataframes.append(df)

# # ---- plot

# font
font_prop = fm.FontProperties(fname='../../fonts/arial.ttf')
plt.rcParams.update({
    'font.family': font_prop.get_name(),
    'font.size': 14})

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
    ntypes = len(df)
    colors = cmap((np.arange(ntypes)/ntypes) * .8)

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

plt.tight_layout()
plt.savefig('pie_chart_for_all.png', dpi=300)
