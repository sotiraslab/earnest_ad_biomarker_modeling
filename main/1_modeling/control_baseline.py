#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:22:26 2023

@author: earnestt1234
"""

# ---- imports

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
OUTPUT = 'control_baseline'

# ---- Read

df = pd.read_csv(PATH_DATA)
df['Sex'] = np.where(df['Sex'] == 'Male', 1., 0.)
df['HasE4'] = df['HasE4'].astype(float)

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

def my_hstack(x):
    return np.hstack(x) if x else x

def get_covariates(fitted_models, data):

    X = []
    for m in fitted_models:
        X.append(m.covariates(data))

    return np.hstack(X) if X else None

def get_training_x(data, models):
    base_x = data[COVARIATES].to_numpy()
    if len(models) == 0:
        return base_x
    else:
        return np.hstack([base_x, get_covariates(models, data)])

def test_atn_model(models, train_data, test_data):

    if not isinstance(models, list):
        models = [models]

    for m in models:
        m.fit(train_data)

    X_train = get_training_x(train_data, models)
    y_train = train_data[TARGET].to_numpy()
    omit = np.any(np.isnan(X_train), axis=1)
    X_train = X_train[~omit, :]
    y_train = y_train[~omit]


    X_test = get_training_x(test_data, models)
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

def get_combo_model(selected_models, amyloid=None, tau=None,
                    neurodegeneration=None):

    df = selected_models

    measure_types = [amyloid, tau, neurodegeneration]
    modalities = ['amyloid', 'tau', 'neurodegeneration']
    output = []

    for measure_type, modality in zip(measure_types, modalities):
        if measure_type is None: continue;
        key = df.loc[(df['atn'] == modality) & (df['measure_type'] == measure_type)]['name'].iloc[0]
        output.append(SEPARATE_ATN_MODELS[modality][measure_type][key])

    return output

all_results = []
cn_results = []
dem_results = []

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
            'Baseline': get_combo_model(selected_models, None, None, None),
            'Binary A': get_combo_model(selected_models, 'binary', None, None),
            'Binary T': get_combo_model(selected_models, None, 'binary', None),
            'Binary N': get_combo_model(selected_models, None, None, 'binary'),
            'All binary': get_combo_model(selected_models, 'binary', 'binary', 'binary'),
            'Categorical A': get_combo_model(selected_models, 'categorical', None, None),
            'Categorical T': get_combo_model(selected_models, None, 'categorical', None),
            'Categorical N': get_combo_model(selected_models, None, None, 'categorical'),
            'All categorical': get_combo_model(selected_models, 'categorical', 'categorical', 'categorical'),
            'Continuous A': get_combo_model(selected_models, 'continuous', None, None),
            'Continuous T': get_combo_model(selected_models, None, 'continuous', None),
            'Continuous N': get_combo_model(selected_models, None, None, 'continuous'),
            'All continuous': get_combo_model(selected_models, 'continuous', 'continuous', 'continuous'),
            }

        all_outer = outer_test
        cn_outer = outer_test.loc[outer_test['Dementia'] == 'No']
        dem_outer = outer_test.loc[outer_test['Dementia'] == 'Yes']


        for data, holder in zip([all_outer, cn_outer, dem_outer], [all_results, cn_results, dem_results]):

            for name, model in final_models.items():
                rmse = test_atn_model(model, outer_train, data)
                row = {'model': name,
                       'rmse': rmse,
                       'fold': i,
                       'repeat': r}
                holder.append(row)

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

def compute_stats(results, baseline = 'Baseline'):
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

def results_boxplot(results, baseline='Baseline', save=None, stats=True, title=None):

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
               ['#FFDDAA'] * 3 +
               ['#F7A934'] +
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
                    fontsize=20, fontweight='bold', color='darkorange')


    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

results_boxplot(all_results, save=os.path.join(OUTPUT, 'all_test_set.png'), title='Accuracy (entire test set)')
results_boxplot(cn_results, save=os.path.join(OUTPUT, 'cn_test_set.png'), title='Accuracy (CDR=0)')
results_boxplot(dem_results, save=os.path.join(OUTPUT, 'dementia_test_set.png'), title='Accuracy (CDR>=0.5)')

#%%

a, b = compute_stats(all_results)
