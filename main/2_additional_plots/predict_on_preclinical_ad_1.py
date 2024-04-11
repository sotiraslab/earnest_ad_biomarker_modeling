#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:15:21 2024

@author: earnestt1234
"""

import os
import pickle
import warnings

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import mean_squared_error

import sys
sys.path.append('../1_modeling')

from helpers import results_boxplot, get_training_x

# ---- load modeling results
def load_results_for_folder(folder, result, check_short=True):

    path = os.path.join(folder, result)
    if not os.path.exists(path) and not check_short:
        raise RuntimeError(f'"{result}" for "{folder}" are missing.')

    if not os.path.exists(path):
        short_folder = folder + '_short'
        short_path = os.path.join(short_folder, result)

        if os.path.exists(short_path):
            warnings.warn(f'The full results are missing for "{folder}", but the short results are present.  Using the latter!',
                          RuntimeWarning)
            path = short_path
        else:
            raise RuntimeError(f'"{result}" for "{folder}" are missing.')

    with open(path, 'rb') as f:
        output = pickle.load(f)

    return output

this_dir = os.path.dirname(os.path.abspath(__file__))
exp2_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp2_combo_atn_models_global_cognition'))
exp1_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp1_svms_global_cognition'))

svms = load_results_for_folder(exp1_folder, 'models.pickle')
atn_predictors = load_results_for_folder(exp2_folder, 'models.pickle')
atn_lms = load_results_for_folder(exp2_folder, 'lm.pickle')

# ---- load dataset

dataset = pd.read_csv('../../outputs/maindata/maindata.csv')

# ---- load CV splitter from experiments

covariates=['Age', 'SexBinary', 'HasE4Binary']
stratify='CDRBinned'
repeats=10
outer_splits=10
outer_seed=0

# outer CV loop
rows = []
for r in range(repeats):

    outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)
    for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(dataset, dataset[stratify])):

        index = (r * outer_splits) + i
        outer_test = dataset.iloc[outer_test_index, :]
        outer_test_cn = outer_test.loc[outer_test['CDRBinned'].eq('0.0'), :]
        gt = outer_test_cn['PHC_GLOBAL']

        # predict with LMs
        for name in atn_predictors.keys():
            atn_model = atn_predictors[name][index]
            lm = atn_lms[name][index]
            X = get_training_x(outer_test_cn, covariates, atn_model)
            preds = lm.predict(X)

            row = {'name': name,
                   'kind': 'lm',
                   'repeat': r,
                   'fold': i,
                   'rmse': mean_squared_error(gt, preds, squared=False)}
            rows.append(row)

        # predict with SVMs
        for name, svm_list in svms.items():
            svm = svm_list[i]
            preds = svm.predict(outer_test_cn)

            row = {'name': name,
                   'kind': 'svm',
                   'repeat': r,
                   'fold': i,
                   'rmse': mean_squared_error(gt, preds, squared=False)}
            rows.append(row)
#%%
results = pd.DataFrame(rows)
results_filter = results[~results['name'].str.contains('PVC')]

palette = (['gray'] +
    ['#FFDDAA'] * 3 +
    ['#F7A934'] +
    ['#E59C9C'] * 3 +
    ['#ba635d'] +
    ['#AAC3E9'] * 3 +
    ['#7DA1D8'] +
    ['#B3E2AD'] * 3 +
    ['#99C494'])

fig, _ = results_boxplot(results_filter, 'name', baseline='Baseline', stats_vs_baseline=True,
                         n_test=len(outer_test_index), n_train=len(outer_train_index),
                         palette=palette)

outfolder = '../../outputs/additional_plots'
if not os.path.isdir(outfolder):
    os.mkdir(outfolder)

fig.savefig(os.path.join(outfolder, 'preclinical_ad_combo_vs_baseline.svg'))
