#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import atn_modeling.atn_predictor_classes
from common import load_results, set_font_properties

# ---- load saved linear models

linear_models = load_results('exp2_combo_atn_models_global_cognition', 'lm.pickle')

all_binary = linear_models['All binary']
all_continuous = linear_models['All continuous']

# hard to work with the categorical, since the number
# of coefficients is equal to the number of categories
# all_categorical = linear_models['All categorical']

# --- make dataframes with coefficients for each model

# note: the order is hardcoded right now
# sklearn doesn't save the feature names for coefficients
# this could be explictly saved in the exp2

coef_order = ['amyloid', 'tau', 'neurodegeneration', 'age', 'sex', 'hasE4']

binary_coefs = pd.DataFrame([np.abs(m['lm'].coef_) for m in all_binary], columns=coef_order)
continuous_coefs = pd.DataFrame([np.abs(m['lm'].coef_) for m in all_continuous], columns=coef_order)

# set font
set_font_properties('arial.ttf')

# plot
def jitter(ys, xcenter, spread):
    xs = np.random.uniform(0, spread/2, size=len(ys))
    half = int(len(ys)/2)
    xs[np.arange(len(xs)) < half] *= -1
    np.random.shuffle(xs)
    xs += xcenter
    return xs

fig, ax = plt.subplots()
alpha = .7

y = binary_coefs['amyloid']
x = jitter(y, 0, .75)
ax.scatter(x, y, marker='o', color='#882255', label='Amyloid', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3, zorder=3)

y = binary_coefs['tau']
x = jitter(y, 1, .75)
ax.scatter(x, y, marker='o', color='#117733', label='Tau', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3)

y = binary_coefs['neurodegeneration']
x = jitter(y, 2, .75)
ax.scatter(x, y, marker='o', color='#332288', label='Neurodegeneration', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3)

y = continuous_coefs['amyloid']
x = jitter(y, 4, .75)
ax.scatter(x, y, marker='o', color='#882255', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3)

y = continuous_coefs['tau']
x = jitter(y, 5, .75)
ax.scatter(x, y, marker='o', color='#117733', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3)

y = continuous_coefs['neurodegeneration']
x = jitter(y, 6, .75)
ax.scatter(x, y, marker='o', color='#332288', alpha=alpha, zorder=2, edgecolor='none')
ax.plot([x.min(), x.max()], [y.mean(), y.mean()], color='black', lw=3)

# plot formatting
ax.set_xticks([1, 5])
ax.set_xticklabels(['Binary', 'Continuous'])
ax.set_ylabel('Coefficients')
ax.legend()

outfile = os.path.join('figures/combo_atn_models_linear_model_coefficients.svg')

plt.tight_layout()
plt.savefig(outfile)
