#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os
import pickle

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import atn_modeling.atn_predictor_classes
from common import load_results

# ---- load saved models from binary combination experiment

linear_models = load_results('exp2_combo_atn_models_global_cognition', 'models.pickle')

all_binary = linear_models['All binary']

# ---- collect values into dataframe

df = pd.DataFrame(
    {'iteration': range(len(all_binary)),
      'amyloid_cutoff': [m[0].cutoff for m in all_binary],
      'amyloid_model': [m[0].nickname for m in all_binary],
      'tau_cutoff': [m[1].cutoff for m in all_binary],
      'tau_model': [m[1].nickname for m in all_binary],})

# must covert centiloid to SUVR
# CL=188.22 × SUMMARY_SUVR − 189.16
def centiloid_to_suvr(centiloid):
    return (centiloid + 189.16) / 188.22

df['amyloid_cutoff_converted'] = np.where(df['amyloid_model'].str.contains('Centiloid'),
                                          centiloid_to_suvr(df['amyloid_cutoff']),
                                          df['amyloid_cutoff'])

# ----- plot

def draw_meanline(ax, mean, y=.8, color='black', offset=0.01, text=None, **kwargs):
    ax.axvline(mean, color=color, label=text, **kwargs)
    # transform = ax.get_xaxis_transform()
    # ax.text(x=mean + .01,
    #         y=y,
    #         s=text,
    #         color=color,
    #         transform=transform)

# base plot
try:
    font_prop = fm.FontProperties(fname='arial.ttf')
    plt.rcParams.update({
        'font.family': font_prop.get_name()})
except:
    pass

fig, (amy_ax, tau_ax) = plt.subplots(nrows=2, sharex=True,
                                      figsize=(14, 6))
plt.rcParams.update({'font.size': 20})
amy_ax.set_title('Amyloid')
tau_ax.set_title('Tau')

# amyloid
ax = amy_ax
x = df['amyloid_cutoff_converted']
color = '#882255'

sns.kdeplot(x ,
            ax=ax,
            color=color,
            fill=True)

mean = x.mean().round(3)
ax.axvline(mean, color=color)
draw_meanline(ax, mean, 0.8, color=color, text=f'Cross-validated ({mean})', lw=4)
draw_meanline(ax, 1.11, 0.6, color='k', text='Landau (1.11)')
draw_meanline(ax, 1.24, 0.6, color='brown', text='Su (1.25)')
draw_meanline(ax, 1.42, 0.6, color='orange', text='Jack RW')
draw_meanline(ax, 1.30, 0.6, color='teal', text='Jack Spec.')
draw_meanline(ax, centiloid_to_suvr(15), 0.6, color='red', text='Centiloid>15')
draw_meanline(ax, centiloid_to_suvr(20), 0.6, color='blue', text='Centiloid>20')
draw_meanline(ax, centiloid_to_suvr(25), 0.6, color='green', text='Centiloid>25')
ax.legend(fontsize=14,
          loc='upper left',
          bbox_to_anchor=(1, 1),
          ncol=2)

#tau
ax = tau_ax
x = df['tau_cutoff']
color = '#117733'

transform = ax.get_xaxis_transform()

sns.kdeplot(x ,
            ax=ax,
            color=color,
            fill=True)
mean = x.mean().round(3)
draw_meanline(ax, mean, 0.8, color=color, text=f'Cross-validated ({mean})', lw=4)
draw_meanline(ax, 1.20, 0.8, color='#de1c1a', text='Jack Sens. (1.20)')
draw_meanline(ax, 1.21, 0.6, color='#329f2b', text='Jack Spec. (1.21)')
draw_meanline(ax, 1.23, 0.4, color='#fe7f06', text='Jack Acc-Young (1.23)')
draw_meanline(ax, 1.33, 0.2, color='#6b409f', text='Jack Acc-Matched. (1.33)')
ax.legend(fontsize=14,
          loc='upper left',
          bbox_to_anchor=(1, 1))

# formatting
tau_ax.set_xlabel('Cutoff (SUVR)')

# save
plt.tight_layout()
fig.savefig('figures/binary_pet_cutoffs.svg')
