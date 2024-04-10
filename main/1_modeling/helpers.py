#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 11:34:10 2023

@author: tom.earnest
"""

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pingouin import ttest
from scipy.stats import t
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

from atn_predictor_instances import get_models_by_nickname

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
                        biomarker_col='biomarker',
                        vartype_col='variable_type',
                        name_col='name'):

    df = selected_models

    variable_types = [amyloid, tau, neurodegeneration]
    biomarkers = ['amyloid', 'taupvc' if taupvc else 'tau', 'neurodegeneration']
    output = []

    for variable_type, biomarker in zip(variable_types, biomarkers ):
        if variable_type is None: continue;
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
    
def test_atn_linear_model(models, covariates, target, train_data, test_data):

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

    scaler = StandardScaler()
    X_train_scale = scaler.fit_transform(X_train)
    X_test_scale = scaler.transform(X_test)

    lm = LinearRegression()
    lm.fit(X_train_scale, y_train)

    y_pred = lm.predict(X_test_scale)
    metrics = {'rmse': mean_squared_error(y_test, y_pred, squared=False),
               'r2': r2_score(y_test, y_pred)}

    return metrics, lm

# ---- stats

def compute_stats_vs_baseline(results, n_train, n_test, baseline = 'Baseline',
                              error_measure=True):
    models = [m for m in results['model'].unique() if m != baseline]

    a = results.loc[results['model'] == baseline, 'rmse']
    t_df = []
    corrected_t_df = []
    for m in models:
        alternative = 'greater' if error_measure else 'less'
        side = 'right' if error_measure else 'left'
        b = results.loc[results['model'] == m, 'rmse']
        t_df.append(ttest(a, b, alternative=alternative, paired=True))
        corrected_t_df.append(nadeau_bengio_test(a, b, n_train, n_test, side=side))

    t_df = pd.concat(t_df)
    t_df['model'] = models
    t_df['p-val-fdr'] = multipletests(t_df['p-val'], method='fdr_bh')[1]

    corrected_t_df = pd.DataFrame(corrected_t_df)
    corrected_t_df['model'] = models
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

    return {'t': tstat, 'p': p, 'dof': dof}

# ---- Plot

def results_boxplot(results, groupby, baseline='Baseline', save=None, stats=True,
                    nadeau_bengio=True, title=None, palette=None,
                    n_train=None, n_test=None, order=None,
                    pivot_index=['fold', 'repeat'], pivot_values='rmse',
                    error_measure=True):

    # data
    if order is None:
        order = results[groupby].unique()
    boxplotdata = results.pivot(index=pivot_index, columns=groupby, values=pivot_values)
    boxplotdata = boxplotdata[order]

    # base plot
    try:
        font_prop = fm.FontProperties(fname='../../fonts/arial.ttf')
        plt.rcParams.update({
            'font.family': font_prop.get_name()})
    except:
        pass
    
    fig, ax = plt.subplots(figsize=(6, 8))
    positions = list(range(len(order)))
    positions.reverse()
    bplot = ax.boxplot(boxplotdata, vert=False, patch_artist=True, positions=positions,
                       sym='o', flierprops={'markerfacecolor':'gray', 'markeredgecolor':'gray'})

    # colors
    if palette is None:
        palette = ['Gray'] * len(order)
    
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
    if baseline:
        baseline_median = np.median(boxplotdata[baseline])
        ax.axvline(baseline_median, color='black', linestyle='dashed', zorder=3)

    # stats
    stats_table = None
    
    if stats and baseline is None:
        raise ValueError('Must select baseline if requesting stats.')
    
    if stats:
        idx = bool(nadeau_bengio)
        stats_table = compute_stats_vs_baseline(results, baseline=baseline, n_train=n_train, n_test=n_test,
                                                error_measure=error_measure)[idx]

        xmin, xmax = ax.get_xlim()
        xrng = xmax - xmin

        for i in range(len(stats_table)):

            p = stats_table.iloc[i, :]['p-val-fdr']
            if p > 0.05:
                continue

            stars = '*'
            stars = '**' if p <= 0.01 else stars
            stars = '***' if p <= 0.001 else stars

            model = stats_table.iloc[i, :][groupby]
            rng = (boxplotdata[model].max() - boxplotdata[model].min())
            x =  boxplotdata[model].max() + 0.12 * rng
            y = len(order) - i - 2
            ax.text(x, y, s=stars, rotation=90, ha='center', va='center',
                    fontsize=20, fontweight='bold', color='darkorange')
            while x >= xmax:
                xmax = xmax + (0.05*xrng)
                ax.set_xlim(xmin, xmax)
            
    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)
        
    return ax.get_figure(), stats_table

def results_boxplot_pairwise(results, groupby, baseline='Baseline',
                             pairs=None, save=None,
                             nadeau_bengio=True, title=None, palette=None,
                             n_train=None, n_test=None, order=None,
                             pivot_index=['fold', 'repeat'], pivot_values='rmse'):

    # data
    if order is None:
        order = results[groupby].unique()
    boxplotdata = results.pivot(index=pivot_index, columns=groupby, values=pivot_values)
    boxplotdata = boxplotdata[order]

    # base plot
    try:
        font_prop = fm.FontProperties(fname='../../fonts/arial.ttf')
        plt.rcParams.update({
            'font.family': font_prop.get_name()})
    except:
        pass
    
    fig, ax = plt.subplots(figsize=(6, 8))
    positions = list(range(len(order)))
    positions.reverse()
    bplot = ax.boxplot(boxplotdata, vert=False, patch_artist=True, positions=positions,
                        sym='o', flierprops={'markerfacecolor':'gray', 'markeredgecolor':'gray'})

    # colors
    if palette is None:
        palette = ['Gray'] * len(order)
    
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
    if baseline:
        baseline_median = np.median(boxplotdata[baseline])
        ax.axvline(baseline_median, color='black', linestyle='dashed', zorder=3)

    # stats
    stats_table = None

    if pairs:
        idx = bool(nadeau_bengio)
        stats_table = compute_stats_pairwise(results, pairs=pairs, n_train=n_train, n_test=n_test, name_col='name')[idx]
        
        xmin, xmax = ax.get_xlim()
        xrng = xmax - xmin

        for i in range(len(stats_table)):

            p = stats_table.iloc[i, :]['p-val-fdr']
            if p > 0.05:
                continue

            stars = '*'
            stars = '**' if p <= 0.01 else stars
            stars = '***' if p <= 0.001 else stars

            nameA = stats_table.loc[i, 'a']
            nameB = stats_table.loc[i, 'b']
            mini =  min(boxplotdata[nameA].min(), boxplotdata[nameB].min())
            maxi = max(boxplotdata[nameA].max(), boxplotdata[nameB].max())
            yA = len(order) - list(order).index(nameA) - 1
            yB = len(order) - list(order).index(nameB) - 1
            rng = (maxi - mini)
            xbar =  maxi + 0.12 * rng
            xtext = xbar + (xrng/15)
            ytext = (yA + yB)/2
            
            ax.plot([xbar, xbar], [yA, yB], color='#35abab')
            ax.text(xtext, ytext, s=stars, rotation=90, ha='center', va='center',
                    fontsize=20, fontweight='bold', color='#35abab')
            while xtext >= xmax:
                xmax = xmax + (0.05*xrng)
                ax.set_xlim(xmin, xmax)
            
    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)
        
    return ax.get_figure(), stats_table