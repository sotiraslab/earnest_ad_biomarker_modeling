#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
# from sklearn.metrics import root_mean_squared_error

from atn_modeling import atn_predictor_classes
from atn_modeling.atn_predictor_instances import ATN_PREDICTORS
from atn_modeling.atn_predictor_classes import MultivariateSVR
from atn_modeling.helpers import test_atn_linear_model_return_predictions
from common import load_results

# ---- load data
df = pd.read_csv('datasets/maindata.csv')

#%% --- get Linear model predictions from models based on experiment

lms = load_results('exp2_comboVSsolo_MMSE', 'models.pickle')
lm_keys = [k for k in lms.keys() if 'PVC' not in k and k != 'Baseline']

# repeat CV routine to get test predictions
outer_seed = 0
outer_splits = 10
repeats = 10
stratify='CDRBinned'
covariates=['Age', 'SexBinary', 'HasE4Binary']
target='DeltaMMSE'

indvl = df.loc[:, ['RID', 'TauID', 'CDRBinned', 'DeltaMMSE', 'DeltaMMSEBootSE']].copy()

for key in lm_keys:
    print(key)
    test_predictions = np.zeros(shape=(len(df), repeats * outer_splits))
    test_predictions[:] = np.nan

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

    average_test_predictions = np.nanmean(test_predictions, axis=1)
    indvl[f'{key}_prediction'] = average_test_predictions
    indvl[f'{key}_residual'] = indvl[f'{key}_prediction'] - indvl['DeltaMMSE']

#%% --- Same for SVMs

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
            # rmse.append(root_mean_squared_error(y_true, y_pred))

    average_test_predictions = np.nanmean(test_predictions, axis=1)
    indvl[f'{key}_prediction'] = average_test_predictions
    indvl[f'{key}_residual'] = indvl[f'{key}_prediction'] - indvl['DeltaMMSE']


#%% Save

outputfolder = 'outputs'
outpath = os.path.join(outputfolder, 'individual_predictions.csv')
indvl.to_csv(outpath, index=False)
