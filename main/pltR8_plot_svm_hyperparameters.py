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
import pandas as pd

from atn_modeling.helpers import results_boxplot
from common import load_results, set_labels_baseline_exp, set_labels_binary_exp

svm = load_results('exp1_svm_MMSE', 'results.csv')
svm = svm[~svm['model'].str.contains('PVC')].copy()

# Note that there is a different search space for C with a linear kernel
# But only the RBF kernel was selected
c_search_space = list(2. ** np.arange(-5, 17, 2))
c_labels = [r"$2^{{{}}}$".format(p) for p in np.arange(-5, 17, 2)]
gamma_search_space = list(2. ** np.arange(-15, 5, 2))
gamma_labels = [r"$2^{{{}}}$".format(p) for p in np.arange(-15, 5, 2)]

models = svm['model'].unique()

for model in models:

    # Plot C
    df = svm[svm['model'].eq(model).copy()]
    C = df.groupby('C')['model'].count()
    C = C.reindex(c_search_space).fillna(0).reset_index()
    x = c_labels
    y = C['model'].values

    fig = plt.figure(figsize=(6, 2))

    plt.bar(x, y, zorder=2, edgecolor='k')
    plt.ylim(0, 100)
    plt.grid(zorder=1)
    plt.ylabel('Count')
    plt.xlabel('C')

    # Plot gamma
    df = svm[svm['model'].eq(model).copy()]
    gamma = df.groupby('gamma')['model'].count()
    gamma = gamma.reindex(gamma_search_space).fillna(0).reset_index()
    x = gamma_labels
    y = gamma['model'].values

    fig = plt.figure(figsize=(6, 2))

    plt.bar(x, y, zorder=2, edgecolor='k', color='firebrick')
    plt.ylim(0, 100)
    plt.grid(zorder=1)
    plt.ylabel('Count')
    plt.xlabel('Gamma')
