#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:34:10 2023

@author: tom.earnest
"""

# was running into weird issue where importing
# sklearn stuff later caused massive slow down
# unable to find good documentation
from sklearn.linear_model import LinearRegression
from sklearn.metrics import root_mean_squared_error, r2_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pingouin import ttest
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

from atn_modeling.atn_predictor_instances import get_models_by_nickname

# ---- modeling helpers

def get_atn_covariates(fitted_models, data):

    X = []
    for m in fitted_models:
        X.append(m.covariates(data))

    return np.hstack(X) if X else None

def get_combo_atn_model(selected_models, atn_predictors_dict,
                        amyloid=None, tau=None,
                        neurodegeneration=None,
                        taupvc=False,
                        csf=False,
                        biomarker_col='biomarker',
                        vartype_col='variable_type',
                        name_col='name'):

    df = selected_models

    variable_types = [amyloid, tau, neurodegeneration]
    if csf:
        biomarkers = ['csf_amyloid', 'csf_tau', 'csf_neurodegeneration']
    else:
        biomarkers = ['amyloid', 'taupvc' if taupvc else 'tau', 'neurodegeneration']
    output = []

    for variable_type, biomarker in zip(variable_types, biomarkers ):
        if variable_type is None:
            continue
        key = df.loc[(df[biomarker_col] == biomarker) & (df[vartype_col] == variable_type)][name_col].iloc[0]
        output += get_models_by_nickname(key)

    return output

def get_training_x(data, covariates, models):
    base_x = data[covariates].to_numpy()
    if len(models) == 0:
        return base_x
    else:
        return np.hstack([base_x, get_atn_covariates(models, data)])

def svm_best_param_lookup(results_table, svm_name, svm_params):

    lookup = results_table.to_dict(orient='records')
    d = {}

    for entry in lookup:
        if entry['name'] == svm_name:
            d = entry
            break

    if not d:
        raise ValueError(f"Cannot find params for model '{svm_name}'")

    params = {k: v for k, v in d.items() if k in svm_params}
    return params

def test_atn_linear_model(models, covariates, target, train_data, test_data,
                          scale=False):

    if not isinstance(models, list):
        models = [models]

    for m in models:
        m.fit(train_data)

    X_train = get_training_x(train_data, covariates, models)
    
    y_train = train_data[target].to_numpy()
    omit = np.any(np.isnan(X_train), axis=1)
    X_train = X_train[~omit, :]
    y_train = y_train[~omit]

    X_test = get_training_x(test_data, covariates, models)
    y_test = test_data[target].to_numpy()
    omit = np.any(np.isnan(X_test), axis=1)
    X_test = X_test[~omit, :]
    y_test = y_test[~omit]

    if scale:
        pipeline = Pipeline([('scaler', StandardScaler()),
                             ('lm', LinearRegression())])
    else:
        pipeline = Pipeline([('lm', LinearRegression())])

    pipeline.fit(X_train.copy(), y_train.copy())
    y_pred = pipeline.predict(X_test.copy())

    metrics = {'rmse': root_mean_squared_error(y_test, y_pred),
               'r2': r2_score(y_test, y_pred)}

    return metrics, pipeline

# ---- stats

def compute_stats_vs_baseline(results, n_train, n_test, baseline = 'Baseline',
                              error_measure=True, name_col='model', value_col='rmse'):
    models = [m for m in results[name_col].unique() if m != baseline]

    a = results.loc[results[name_col] == baseline, value_col]
    t_df = []
    corrected_t_df = []
    for m in models:
        alternative = 'greater' if error_measure else 'less'
        side = 'right' if error_measure else 'left'
        b = results.loc[results[name_col] == m, value_col]
        t_df.append(ttest(a, b, alternative=alternative, paired=True))
        corrected_t_df.append(nadeau_bengio_test(a, b, n_train, n_test, side=side))

    t_df = pd.concat(t_df)
    t_df[name_col] = models
    t_df['p-val-fdr'] = multipletests(t_df['p-val'], method='fdr_bh')[1]

    corrected_t_df = pd.DataFrame(corrected_t_df)
    corrected_t_df[name_col] = models
    corrected_t_df['p-val-fdr'] = multipletests(corrected_t_df['p'], method='fdr_bh')[1]

    return t_df, corrected_t_df

def compute_stats_pairwise(results, pairs, n_train, n_test, value='rmse',
                           name_col='model'):

    t_df = []
    corrected_t_df = []
    a_list = [pair[0] for pair in pairs]
    b_list = [pair[1] for pair in pairs]
    for colA, colB in pairs:
        a = results.loc[results[name_col] == colA, value]
        b = results.loc[results[name_col] == colB, value]
        t_df.append(ttest(a, b, alternative='two-sided', paired=True))
        corrected_t_df.append(nadeau_bengio_test(a, b, n_train, n_test, side='both'))

    t_df = pd.concat(t_df)
    t_df['a'] = a_list
    t_df['b'] = b_list
    t_df['p-val-fdr'] = multipletests(t_df['p-val'], method='fdr_bh')[1]

    corrected_t_df = pd.DataFrame(corrected_t_df)
    corrected_t_df['a'] = a_list
    corrected_t_df['b'] = b_list
    corrected_t_df['p-val-fdr'] = multipletests(corrected_t_df['p'], method='fdr_bh')[1]

    return t_df, corrected_t_df


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

    return {'mean_a': np.mean(a),
            'sd_a': np.std(a),
            'mean_b': np.mean(b),
            'sd_b': np.std(b),
            'mean_diff': numerator,
            't': tstat,
            'p': p,
            'dof': dof}

# ---- Plot

def results_boxplot(results, groupby, baseline='Baseline', save=None,
                    stats_vs_baseline=False, stats_pairs=None,
                    stats_pairs_positions=None,
                    nadeau_bengio=True, title=None, palette=None,
                    n_train=None, n_test=None, order=None,
                    pivot_index=['fold', 'repeat'], pivot_values='rmse',
                    error_measure=True, hatch=None, font_file=None):

    # data
    if order is None:
        order = results[groupby].unique()
    boxplotdata = results.pivot(index=pivot_index, columns=groupby, values=pivot_values)
    boxplotdata = boxplotdata[order]

    # font
    if font_file:
        font_prop = fm.FontProperties(fname=font_file)
        val = font_prop.get_name()
    else:
        val = 'Arial'
    plt.rcParams.update({'font.family': val,
                         'mathtext.fontset': 'custom',
                         'mathtext.rm': val,
                         'mathtext.it': f'{val}:italic',
                         'mathtext.bf': f'{val}:bold',
                         'mathtext.cal': val
                         })

    # base plot
    fig, ax = plt.subplots(figsize=(6, 8))
    positions = list(range(len(order)))
    positions.reverse()
    bplot = ax.boxplot(boxplotdata, vert=False, patch_artist=True, positions=positions,
                       sym='o', flierprops={'markerfacecolor':'gray', 'markeredgecolor':'gray'})

    # colors
    if palette is None:
        palette = ['Gray'] * len(order)

    if hatch is None:
        hatch = [False] * len(order)

    double_palette = np.repeat(palette, 2)

    for patch, color, h in zip(bplot['boxes'], palette, hatch):
        if h:
            patch.set_facecolor(color)
            patch.set_edgecolor('white')
            patch.set_hatch('//')
        else:
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
    if baseline:
        baseline_median = np.median(boxplotdata[baseline])
        ax.axvline(baseline_median, color='black', linestyle='dashed', zorder=3)

    # stats
    stats_tbl_bl = None
    stats_tbl_pairs = None

    if stats_vs_baseline and baseline is None:
        raise ValueError('Must select baseline if requesting stats.')

    if stats_vs_baseline:
        idx = bool(nadeau_bengio)
        stats_tbl_bl = compute_stats_vs_baseline(results, baseline=baseline, n_train=n_train, n_test=n_test,
                                                 error_measure=error_measure,
                                                 name_col=groupby)[idx]

        xmin, xmax = ax.get_xlim()
        xrng = xmax - xmin

        for i in range(len(stats_tbl_bl)):

            p = stats_tbl_bl.iloc[i, :]['p-val-fdr']
            if p > 0.05:
                continue

            stars = '*'
            stars = '**' if p <= 0.01 else stars
            stars = '***' if p <= 0.001 else stars

            model = stats_tbl_bl.iloc[i, :][groupby]
            rng = (boxplotdata[model].max() - boxplotdata[model].min())
            x =  boxplotdata[model].max() + 0.12 * rng
            y = len(order) - list(order).index(model) - 1
            ax.text(x, y, s=stars, rotation=90, ha='center', va='center',
                    fontsize=16, fontweight='bold', color='darkorange')
            while x >= xmax:
                xmax = xmax + (0.05*xrng)
                ax.set_xlim(xmin, xmax)

    if stats_pairs:
        idx = bool(nadeau_bengio)
        stats_tbl_pairs = compute_stats_pairwise(results,
                                                 pairs=stats_pairs,
                                                 n_train=n_train,
                                                 n_test=n_test,
                                                 name_col=groupby)[idx]

        xmin, xmax = ax.get_xlim()
        xrng = xmax - xmin

        for i in range(len(stats_tbl_pairs)):

            p = stats_tbl_pairs.iloc[i, :]['p-val-fdr']
            if p > 0.05:
                continue

            stars = '*'
            stars = '**' if p <= 0.01 else stars
            stars = '***' if p <= 0.001 else stars

            nameA = stats_tbl_pairs.loc[i, 'a']
            nameB = stats_tbl_pairs.loc[i, 'b']
            mini =  min(boxplotdata[nameA].min(), boxplotdata[nameB].min())
            maxi = max(boxplotdata[nameA].max(), boxplotdata[nameB].max())
            yA = len(order) - list(order).index(nameA) - 1
            yB = len(order) - list(order).index(nameB) - 1
            rng = (maxi - mini)
            xbar =  (maxi + 0.15 * rng) if stats_pairs_positions is None else stats_pairs_positions[i]
            xtext = xbar + (xrng/20)
            ytext = (yA + yB)/2

            tip = xrng/100
            clr = 'dimgray'
            ax.plot([xbar, xbar], [yA, yB], color=clr)
            ax.plot([xbar-tip, xbar], [yA, yA], color=clr)
            ax.plot([xbar-tip, xbar], [yB, yB], color=clr)
            ax.text(xtext, ytext, s=stars, rotation=90, ha='center', va='center',
                    fontsize=16, fontweight='bold', color=clr)
            while xtext >= xmax:
                xmax = xmax + (0.05*xrng)
                ax.set_xlim(xmin, xmax)

    stats_dict = {'baseline': stats_tbl_bl,
                  'pairs': stats_tbl_pairs}

    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)

    return ax.get_figure(), stats_dict
