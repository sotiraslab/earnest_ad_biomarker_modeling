#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os

import matplotlib.pyplot as plt
import numpy as np

from common import load_results, set_font_properties

set_font_properties()

svm = load_results('exp1_svm_MMSE', 'results.csv')
svm = svm[~svm['model'].str.contains('PVC')].copy()

# Note that there is a different search space for C with a linear kernel
# But only the RBF kernel was selected
c_search_space = list(2. ** np.arange(-5, 17, 2))
c_labels = [r"$2^{{{}}}$".format(p) for p in np.arange(-5, 17, 2)]
gamma_search_space = list(2. ** np.arange(-15, 5, 2))
gamma_labels = [r"$2^{{{}}}$".format(p) for p in np.arange(-15, 5, 2)]

models = svm['model'].unique()

odir = os.path.join('figures', 'svm_hyperparameters')
os.makedirs(odir, exist_ok=True)

for model in models:

    model_name = model.lower().replace(' ', '_')

    # Plot C
    df = svm[svm['model'].eq(model).copy()]
    C = df.groupby('C')['model'].count()
    C = C.reindex(c_search_space).fillna(0).reset_index()
    x = c_labels
    y = C['model'].values

    fig = plt.figure(figsize=(3.5, 1.5))

    plt.bar(x, y, zorder=2, edgecolor='k')
    plt.ylim(0, 100)
    plt.grid(zorder=1)
    plt.ylabel('Count')
    plt.xlabel('C')

    plotpath = os.path.join(odir, f'{model_name}_C.svg')
    plt.savefig(plotpath)

    # Plot gamma
    df = svm[svm['model'].eq(model).copy()]
    gamma = df.groupby('gamma')['model'].count()
    gamma = gamma.reindex(gamma_search_space).fillna(0).reset_index()
    x = gamma_labels
    y = gamma['model'].values

    fig = plt.figure(figsize=(3.5, 1.5))

    plt.bar(x, y, zorder=2, edgecolor='k', color='firebrick')
    plt.ylim(0, 100)
    plt.grid(zorder=1)
    plt.ylabel('Count')
    plt.xlabel('Gamma')

    plotpath = os.path.join(odir, f'{model_name}_gamma.svg')
    plt.savefig(plotpath)
