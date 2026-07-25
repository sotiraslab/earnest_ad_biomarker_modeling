#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os

from scipy.stats import spearmanr
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
# from sklearn.metrics import root_mean_squared_error

from atn_modeling import atn_predictor_classes
from atn_modeling.atn_predictor_instances import ATN_PREDICTORS
from atn_modeling.atn_predictor_classes import MultivariateSVR
from atn_modeling.helpers import test_atn_linear_model_return_predictions, results_boxplot
from common import load_results

# ---- load data
df = pd.read_csv('datasets/maindata.csv')


#%% --- LMs: combination ATN vs baseline

spearman_combo_v_baseline = {}

lms = load_results('exp2_comboVSsolo_MMSE', 'models.pickle')
lm_keys = [k for k in lms.keys() if 'PVC' not in k]

# repeat CV routine to get test predictions
outer_seed = 0
outer_splits = 10
repeats = 10
stratify='CDRBinned'
covariates=['Age', 'SexBinary', 'HasE4Binary']
target='DeltaMMSE'

for key in lm_keys:
    print(key)
    test_predictions = np.zeros(shape=(len(df), repeats * outer_splits))
    test_predictions[:] = np.nan

    spearman_combo_v_baseline[key] = []

    for r in range(repeats):

        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)

        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(df, df[stratify])):
            iteration = (r * outer_splits) + i
            outer_train = df.iloc[outer_train_index, :]
            outer_test = df.iloc[outer_test_index, :]

            model = lms[key][iteration]
            metrics, y_test, y_pred = test_atn_linear_model_return_predictions(
                models=model,
                covariates=covariates,
                target=target,
                train_data=outer_train,
                test_data=outer_test)
            test_predictions[outer_test.index, iteration] = y_pred
            spearman_combo_v_baseline[key].append(metrics['spearman'])

#%% --- LMs: combination ATN vs binary

spearman_combo_v_binary = {}

lms = load_results('exp3_comboVSbinary_MMSE', 'models.pickle')
lm_keys = [k for k in lms.keys() if 'PVC' not in k]

# repeat CV routine to get test predictions
outer_seed = 0
outer_splits = 10
repeats = 10
stratify='CDRBinned'
covariates=['Age', 'SexBinary', 'HasE4Binary']
target='DeltaMMSE'

for key in lm_keys:
    print(key)
    test_predictions = np.zeros(shape=(len(df), repeats * outer_splits))
    test_predictions[:] = np.nan

    spearman_combo_v_binary[key] = []

    for r in range(repeats):

        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)

        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(df, df[stratify])):
            iteration = (r * outer_splits) + i
            outer_train = df.iloc[outer_train_index, :]
            outer_test = df.iloc[outer_test_index, :]

            model = lms[key][iteration]
            metrics, y_test, y_pred = test_atn_linear_model_return_predictions(
                models=model,
                covariates=covariates,
                target=target,
                train_data=outer_train,
                test_data=outer_test)
            test_predictions[outer_test.index, iteration] = y_pred
            spearman_combo_v_binary[key].append(metrics['spearman'])

#%% --- Same for SVMs

spearman_svm = {}

svms = load_results('exp1_svm_MMSE', 'models.pickle')
svm_keys = [key for key in svms.keys() if 'PVC' not in key]

target='DeltaMMSE'
covariates=['Age', 'SexBinary', 'HasE4Binary']
stratify='CDRBinned'
repeats=10
outer_splits=10
outer_seed=0

test_predictions = np.zeros(shape=(len(df), repeats * outer_splits))
test_predictions[:] = np.nan

key = 'ATM_SVM'

for key in svm_keys:

    spearman_svm[key] = []

    for r in range(repeats):
        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)

        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(df, df[stratify])):
            iteration = (r * outer_splits) + i
            outer_train = df.iloc[outer_train_index, :]
            outer_test = df.iloc[outer_test_index, :]

            model = svms[key][iteration]
            X = outer_test[model.predictors].values
            y_pred = model.pipeline.predict(X)
            y_true = df.loc[outer_test.index, 'DeltaMMSE']

            test_predictions[outer_test.index, iteration] = y_pred
            spearman_svm[key].append(spearmanr(y_true, y_pred).statistic)

#%% Save

odir = os.path.join('figures', 'spearman_correlations')
os.makedirs(odir, exist_ok=True)

df_combo_v_baseline = pd.DataFrame(spearman_combo_v_baseline)
df_combo_v_baseline.to_csv(os.path.join(odir, 'combo_vs_baseline.csv'), index=False)

df_combo_v_binary = pd.DataFrame(spearman_combo_v_binary)
df_combo_v_binary.to_csv(os.path.join(odir, 'combo_vs_binary.csv'), index=False)

df_svm = pd.DataFrame(spearman_svm)
df_svm.to_csv(os.path.join(odir, 'svm.csv'), index=False)
