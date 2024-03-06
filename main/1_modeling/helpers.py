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
from statsmodels.stats.multitest import multipletests

# ---- modeling helpers

def get_atn_covariates(fitted_models, data):

    X = []
    for m in fitted_models:
        X.append(m.covariates(data))

    return np.hstack(X) if X else None

def get_combo_atn_model(selected_models, atn_predictors_dict,
                        amyloid=None, tau=None,
                        neurodegeneration=None):

    df = selected_models

    measure_types = [amyloid, tau, neurodegeneration]
    modalities = ['amyloid', 'tau', 'neurodegeneration']
    output = []

    for measure_type, modality in zip(measure_types, modalities):
        if measure_type is None: continue;
        key = df.loc[(df['atn'] == modality) & (df['measure_type'] == measure_type)]['name'].iloc[0]
        output.append(atn_predictors_dict[modality][measure_type][key])

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

    lm = LinearRegression()
    lm.fit(X_train, y_train)

    y_pred = lm.predict(X_test)
    metrics = {'rmse': mean_squared_error(y_test, y_pred, squared=False),
               'r2': r2_score(y_test, y_pred)}

    return metrics

# ---- stats

def compute_stats(results, n_train, n_test, baseline = 'Baseline'):
    models = [m for m in results['model'].unique() if m != baseline]

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
                    n_train=None, n_test=None, order=None):

    # data
    if order is None:
        order = results[groupby].unique()
    boxplotdata = results.pivot(index=['fold', 'repeat'], columns=groupby, values='rmse')
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
    baseline_median = np.median(boxplotdata[baseline])
    ax.axvline(baseline_median, color='black', linestyle='dashed', zorder=3)

    # stats
    stats_table = None
    if stats:
        idx = bool(nadeau_bengio)
        stats_table = compute_stats(results, baseline=baseline, n_train=n_train, n_test=n_test)[idx]

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


    # save
    if save is not None:
        plt.tight_layout()
        fig.savefig(save, dpi=300)
        
    return ax.get_figure(), stats_table