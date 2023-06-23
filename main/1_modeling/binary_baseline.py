#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:22:26 2023

@author: earnestt1234
"""

# ---- imports

from copy import deepcopy
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pingouin import ttest
from scipy.stats import t
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import StratifiedKFold
from statsmodels.stats.multitest import multipletests


from model_classes import (BinaryAsIs,
                           BinaryGMM,
                           BinaryManual,
                           BinaryZScore,
                           CategoricalStager,
                           Continuous,
                           MultivariateSVR,
                           Quantiles)

# ---- Variables
OUTER_SPLITS = 5
INNER_SPLITS = 5
REPEATS = 20

TARGET = 'PACC.ADNI'
COVARIATES = ['Age', 'Sex', 'HasE4']
STRATIFY_COLUMN = 'CDRBinned'

# ---- Paths

PATH_DATA = '../../data/derivatives/adni_base_table.csv'
OUTPUT = 'binary_baseline'

# ---- Read

df = pd.read_csv(PATH_DATA)
df['Sex'] = np.where(df['Sex'] == 'Male', 1., 0.)
df['HasE4'] = df['HasE4'].astype(float)
# df.loc[df['Dementia'] == 'Unknown', 'Dementia'] = 'No'

df = df[~df['CDRBinned'].isna()]
df = df[~df['PACC.ADNI'].isna()]
df = df[~df['HasE4'].isna()]

amy_columns = list(df.columns[df.columns.str.startswith('AV45')])
tau_columns = list(df.columns[df.columns.str.startswith('FTP')])
gm_columns = list(df.columns[df.columns.str.endswith('VOLUME') & ~ df.columns.str.contains('BRAAK|META')])
roi_columns = amy_columns + tau_columns + gm_columns

amy_columns += COVARIATES
tau_columns += COVARIATES
gm_columns += COVARIATES
roi_columns += COVARIATES

# ---- Setup models

SEPARATE_ATN_MODELS = {
    'amyloid': {
        'binary': {
            'ADNI_analysis': BinaryAsIs('AmyloidPositive'),
            'centiloid_20': BinaryManual('Centiloid', 20)},
        'categorical': {
            'centiloid_quantiles' : Quantiles('Centiloid')},
        'continuous': {
            'centiloid': Continuous('Centiloid')}},
    'tau': {
        'binary': {
            'mtt_gmm': BinaryGMM('META_TEMPORAL_SUVR'),
            'mtt_z2.0': BinaryZScore('META_TEMPORAL_SUVR', 'Control', cutoff=2.0),
            'mtt_z2.5': BinaryZScore('META_TEMPORAL_SUVR', 'Control', cutoff=2.5),
            'mtt_1.20': BinaryManual('META_TEMPORAL_SUVR', cutoff=1.20),
            'mtt_1.21': BinaryManual('META_TEMPORAL_SUVR', cutoff=1.21),
            'mtt_1.23': BinaryManual('META_TEMPORAL_SUVR', cutoff=1.23),
            'mtt_1.33':  BinaryManual('META_TEMPORAL_SUVR', cutoff=1.33)},
        'categorical': {
            'mtt_quantiles': Quantiles('META_TEMPORAL_SUVR'),
            'braak_stage_gmm': CategoricalStager(['BRAAK1_SUVR', 'BRAAK34_SUVR', 'BRAAK56_SUVR'])},
        'continuous': {
            'mtt': Continuous('META_TEMPORAL_SUVR'),
            'braak1': Continuous('BRAAK1_SUVR'),
            'braak34': Continuous('BRAAK34_SUVR'),
            'braak56': Continuous('BRAAK56_SUVR')}
        },
    'neurodegeneration': {
        'binary': {
            'hipp_z2.0': BinaryZScore('HIPPOCAMPUS_VOLUME', control_col='Control', cutoff=-2.0, greater=False),
            'hipp_z2.5': BinaryZScore('HIPPOCAMPUS_VOLUME', control_col='Control', cutoff=-2.5, greater=False),
            'mttvol_z2.0': BinaryZScore('META_TEMPORAL_VOLUME', control_col='Control', cutoff=-2.0, greater=False),
            'mttvol_z2.5': BinaryZScore('META_TEMPORAL_VOLUME', control_col='Control', cutoff=-2.5, greater=False)},
        'categorical': {
            'hipp_quantiles': Quantiles('HIPPOCAMPUS_VOLUME'),
            'mttvol_quantiles': Quantiles('META_TEMPORAL_VOLUME')},
        'continuous': {
            'hipp': Continuous('HIPPOCAMPUS_VOLUME'),
            'mttvol': Continuous('META_TEMPORAL_VOLUME')}},
    }

SVM_MODELS = {'SVM (amyloid)': MultivariateSVR(amy_columns, TARGET),
              'SVM (tau)': MultivariateSVR(tau_columns, TARGET),
              'SVM (GM)': MultivariateSVR(gm_columns, TARGET),
              'SVM (combined)': MultivariateSVR(roi_columns, TARGET)}

# ---- modeling helpers

def get_covariates(fitted_models, data):

    X = []
    for m in fitted_models:
        X.append(m.covariates(data))

    return np.hstack(X)

def test_atn_model(models, train_data, test_data):

    if not isinstance(models, list):
        models = [models]

    for m in models:
        m.fit(train_data)

    X_train = np.hstack([train_data[COVARIATES].to_numpy(), get_covariates(models, train_data)])
    y_train = train_data[TARGET].to_numpy()
    omit = np.any(np.isnan(X_train), axis=1)
    X_train = X_train[~omit, :]
    y_train = y_train[~omit]


    X_test = np.hstack([test_data[COVARIATES].to_numpy(), get_covariates(models, test_data)])
    y_test = test_data[TARGET].to_numpy()
    omit = np.any(np.isnan(X_test), axis=1)
    X_test = X_test[~omit, :]
    y_test = y_test[~omit]

    lm = LinearRegression()
    lm.fit(X_train, y_train)

    y_pred = lm.predict(X_test)
    rmse = mean_squared_error(y_test, y_pred, squared=False)

    return rmse

# ---- Main

def get_combo_model(selected_models, amyloid='binary', tau='binary',
                    neurodegeneration='binary'):

    df = selected_models
    amy_key = df.loc[(df['atn'] == 'amyloid') & (df['measure_type'] == amyloid)]['name'].iloc[0]
    tau_key = df.loc[(df['atn'] == 'tau') & (df['measure_type'] == tau)]['name'].iloc[0]
    ndg_key = df.loc[(df['atn'] == 'neurodegeneration') & (df['measure_type'] == neurodegeneration)]['name'].iloc[0]

    return [SEPARATE_ATN_MODELS['amyloid'][amyloid][amy_key],
            SEPARATE_ATN_MODELS['tau'][tau][tau_key],
            SEPARATE_ATN_MODELS['neurodegeneration'][neurodegeneration][ndg_key]]


all_results = []
cn_results = []
dem_results = []

saved_binary_models = []

for r in range(REPEATS):

    print(f"CV Repeat #{r}...")

    outer_cv = StratifiedKFold(n_splits=OUTER_SPLITS, random_state=r+200, shuffle=True)
    inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, random_state=r+2000, shuffle=True)

    for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(df, df[STRATIFY_COLUMN])):
        outer_train = df.iloc[outer_train_index, :]
        outer_test = df.iloc[outer_test_index, :]

        inner_cv_results = []

        for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[STRATIFY_COLUMN])):
            inner_train = df.iloc[inner_train_index, :]
            inner_test = df.iloc[inner_test_index, :]

            for atn, measure_dict in SEPARATE_ATN_MODELS.items():
                for measure_type, model_dict in measure_dict.items():
                    for name, model in model_dict.items():
                        rmse = test_atn_model(model, inner_train, inner_test)
                        row = {'atn': atn,
                               'measure_type': measure_type,
                               'name': name,
                               'fold': j,
                               'rmse': rmse}
                        inner_cv_results.append(row)

        inner_cv_results = pd.DataFrame(inner_cv_results)
        model_averages = inner_cv_results.groupby(['atn', 'measure_type', 'name'])['rmse'].agg(mean=np.mean, std=np.std).reset_index()
        best_by_measure = model_averages.groupby(['atn', 'measure_type'])['mean'].idxmin()
        selected_models = model_averages.iloc[best_by_measure]

        final_models = {
            'All binary': get_combo_model(selected_models, 'binary', 'binary', 'binary'),
            'Categorical A': get_combo_model(selected_models, 'categorical', 'binary', 'binary'),
            'Categorical T': get_combo_model(selected_models, 'binary', 'categorical', 'binary'),
            'Categorical N': get_combo_model(selected_models, 'binary', 'binary', 'categorical'),
            'All categorical': get_combo_model(selected_models, 'categorical', 'categorical', 'categorical'),
            'Continuous A': get_combo_model(selected_models, 'continuous', 'binary', 'binary'),
            'Continuous T': get_combo_model(selected_models, 'binary', 'continuous', 'binary'),
            'Continuous N': get_combo_model(selected_models, 'binary', 'binary', 'continuous'),
            'All continuous': get_combo_model(selected_models, 'continuous', 'continuous', 'continuous'),
            }

        all_outer = outer_test
        cn_outer = outer_test.loc[outer_test['Dementia'] == 'No']
        dem_outer = outer_test.loc[outer_test['Dementia'] == 'Yes']


        for counter, (data, holder) in enumerate(zip([all_outer, cn_outer, dem_outer], [all_results, cn_results, dem_results])):

            for name, model in final_models.items():
                rmse = test_atn_model(model, outer_train, data)
                row = {'model': name,
                       'rmse': rmse,
                       'fold': i,
                       'repeat': r}
                holder.append(row)
                if name == 'All binary' and counter == 0:
                    saved_binary_models.append(deepcopy(model))

            for name, model in SVM_MODELS.items():
                model.fit(outer_train)
                preds = model.predict(data)
                row = {'model': name,
                       'rmse': mean_squared_error(data[TARGET], preds, squared=False),
                       'fold': i,
                       'repeat': r}
                holder.append(row)

all_results = pd.DataFrame(all_results)
cn_results = pd.DataFrame(cn_results)
dem_results = pd.DataFrame(dem_results)

#%% stats

def nadeau_bengio_test(a, b, n_train, n_test, alpha=0.05, side='both'):
    """
    Implementation of the Nadeau Bengio correction
    for the t-test [1].  This is recommended for accounting
    for the dependence of data points when comparing cross-validated
    model performance.

    Formula is based on the equation outlined by
    Bouckaert & Frank [2]; see section 3.2.

    Implementation also follows this Gist [3], but that
    has some errors that are fixed here.

    [1] https://proceedings.neurips.cc/paper/1999/hash/7d12b66d3df6af8d429c1a357d8b9e1a-Abstract.html
    [2] https://link.springer.com/chapter/10.1007/978-3-540-24775-3_3
    [3] https://gist.github.com/jensdebruijn/13e8eeda85eb8644ac2a4ac4c3b8e732

    Parameters
    ----------
    a : array
        Model A performance metrics.
    b : array
        Model B performance metrics.
    n_train : int
        Number of observations used in a single training fold.
    n_test : int
        Number of observations used in a single testing fold.

    Returns
    -------
    None.

    """

    # check arguments
    if len(a) != len(b):
        raise ValueError("`a` and `b` inputs must be arrays of same length.")

    if side not in ['left', 'right', 'both']:
        raise ValueError("`side` must be 'left', 'right', or 'both'")

    # set variables for equation
    x = np.array(a) - np.array(b)
    var = np.var(x, ddof=1)
    n = len(x)
    n2 = n_test
    n1 = n_train

    # calculate statistic
    numerator = np.mean(x)
    parens = ((1/n) + (n2/n1))
    denominator = np.sqrt(parens * var)
    tstat = numerator / denominator

    # calculate p-value
    dof = n - 1

    if side == 'left':
        p = t.cdf(tstat, dof)
    elif side == 'right':
        p = 1 - t.cdf(tstat, dof)
    elif side == 'both':
        p = 2 * (1 - t.cdf(abs(tstat), dof))

    return {'t': tstat, 'p': p, 'dof': dof}

def compute_stats(results, baseline = 'All binary'):
    models = [m for m in results['model'].unique() if m != baseline]
    n_test = int((1/OUTER_SPLITS) * len(df))
    n_train = len(df) - n_test

    a = results.loc[results['model'] == baseline, 'rmse']
    t_df = []
    corrected_t_df = []
    for m in models:
        b = results.loc[results['model'] == m, 'rmse']
        t_df.append(ttest(a, b, alternative='greater', paired=True))
        corrected_t_df.append(nadeau_bengio_test(a, b, n_train, n_test, side='right'))

    t_df = pd.concat(t_df)
    t_df['model'] = models
    t_df['p-val-fdr'] = multipletests(t_df['p-val'], method='fdr_bh')[1]

    corrected_t_df = pd.DataFrame(corrected_t_df)
    corrected_t_df['model'] = models
    corrected_t_df['p-val-fdr'] = multipletests(corrected_t_df['p'], method='fdr_bh')[1]

    return t_df, corrected_t_df

#%% plot

# style
plt.rcParams.update({'font.family': 'Arial',
                     'font.size': 15})

def results_boxplot(results, baseline='All binary', save=None, stats=True, title=None):

    # data
    order = (list(final_models.keys()) + list(SVM_MODELS.keys()))
    boxplotdata = results.pivot(index=['fold', 'repeat'],columns='model', values='rmse')
    boxplotdata = boxplotdata[order]

    # base plot
    fig, ax = plt.subplots(figsize=(6, 8))
    positions = list(range(len(order)))
    positions.reverse()
    bplot = ax.boxplot(boxplotdata, vert=False, patch_artist=True, positions=positions,
                       sym='o', flierprops={'markerfacecolor':'gray', 'markeredgecolor':'gray'})

    # colors
    palette = (['gray'] +
               ['#E59C9C'] * 3 +
               ['#ba635d'] +
               ['#AAC3E9'] * 3 +
               ['#7DA1D8'] +
               ['#B3E2AD'] * 3 +
               ['#99C494'])

    double_palette = np.repeat(palette, 2)

    for patch, color in zip(bplot['boxes'], palette):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)

    for whiskers, color in zip(bplot['whiskers'], double_palette):
        whiskers.set_color(color)

    for cap, color in zip(bplot['caps'], double_palette):
        cap.set_color(color)

    for median in bplot['medians']:
        median.set_color('white')

    # labels
    ax.set_yticklabels(order)
    ax.set_xlabel('RMSE')
    ax.grid(alpha=.4)

    if title:
        ax.set_title(title, loc='left')

    # baseline
    baseline_median = np.median(boxplotdata[baseline])
    ax.axvline(baseline_median, color='black', linestyle='dashed', zorder=3)

    # stats
    if stats is not None:
        stats = compute_stats(results, baseline=baseline)[1]

        for i in range(len(stats)):

            p = stats.iloc[i, :]['p-val-fdr']
            if p > 0.05:
                continue

            stars = '*'
            stars = '**' if p <= 0.01 else stars
            stars = '***' if p <= 0.001 else stars

            model = stats.iloc[i, :]['model']
            rng = (boxplotdata[model].max() - boxplotdata[model].min())
            x =  boxplotdata[model].max() + 0.12 * rng
            y = len(order) - i - 2
            ax.text(x, y, s=stars, rotation=90, ha='center', va='center',
                    fontsize=25, fontweight='bold', color='darkorange')


    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

results_boxplot(all_results, save=os.path.join(OUTPUT, 'all_test_set.png'), title='Accuracy (entire test set)')
results_boxplot(cn_results, save=os.path.join(OUTPUT, 'cn_test_set.png'), title='Accuracy (CDR=0)')
results_boxplot(dem_results, save=os.path.join(OUTPUT, 'dementia_test_set.png'), title='Accuracy (CDR>=0.5)')

#%% binary selected models

a_models = [m[0] for m in saved_binary_models]
t_models = [m[1] for m in saved_binary_models]
n_models = [m[2] for m in saved_binary_models]

#%% pie chart of selected models

# rename
a_names = []
t_names = []
n_names = []

for m in a_models:
    name = m.__class__.__name__

    if name == 'BinaryAsIs':
        new = 'ADNI assignment'
    elif name == 'BinaryManual':
        new = 'Centiloid >= 20'
    else:
        raise ValueError(name)

    a_names.append(new)

for m in t_models:
    name = m.__class__.__name__

    if name == 'BinaryGMM':
        new = f'GMM({m.y_col})'
    elif name == 'BinaryZScore':
        new = f'{m.cutoff} Z-score ({m.y_col})'
    else:
        raise ValueError(name)

    t_names.append(new)

for m in n_models:
    name = m.__class__.__name__

    if name == 'BinaryZScore':
        new = f'{m.cutoff} Z-score ({m.y_col})'
    else:
        raise ValueError(name)

    n_names.append(new)

from collections import Counter

a_count = Counter(a_names)
t_count = Counter(t_names)
n_count = Counter(n_names)

plt.rcParams.update({'font.family': 'Arial',
                     'font.size': 10})

fig, axes = plt.subplots(1, 3, figsize=(12, 6), dpi=300)

ax = axes[0]
count = a_count
title = 'Amyloid'
wedges, texts = ax.pie(count.values())
ax.legend(wedges, count.keys(), loc='upper center', bbox_to_anchor=(0.5, 0.00))
ax.set_title(title)

ax = axes[1]
count = t_count
title = 'Tau'
wedges, texts = ax.pie(count.values())
ax.legend(wedges, count.keys(), loc='upper center', bbox_to_anchor=(0.5, 0.00))
ax.set_title(title)

ax = axes[2]
count = n_count
title = 'Neurodegeneration'
wedges, texts = ax.pie(count.values())
ax.legend(wedges, count.keys(), loc='upper center', bbox_to_anchor=(0.5, 0.00))
ax.set_title(title)

plt.tight_layout()
fig.savefig('method_pie_chart.png')

#%% tau cutoffs

import seaborn as sns

cutoffs = []

for m in t_models:
    name = m.__class__.__name__

    if name == 'BinaryGMM':
        c = m.cutoff
    elif name == 'BinaryZScore':
        c = (m.cutoff * m.std) + m.mean
    else:
        raise ValueError(name)

    cutoffs.append(c)

fig, ax = plt.subplots(figsize=(8, 3))
sns.histplot(cutoffs, kde=True, ax=ax, label='learned tau cutoffs')
ax.axvline(1.21, color='green', label='Jack (Spec.)')
ax.axvline(1.20, color='firebrick', label='Jack (Sens.)')
ax.axvline(1.23, color='orange', label='Jack (Acc-Young)')
ax.axvline(1.33, color='navy', label='Jack (Acc-Matched)')
ax.set_xlabel('SUVR')

ax.legend(bbox_to_anchor=(1, 1))

plt.tight_layout()
plt.savefig('tau_cutoffs.png', dpi=300)
