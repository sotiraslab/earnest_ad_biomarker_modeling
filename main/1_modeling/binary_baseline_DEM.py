#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:22:26 2023

@author: earnestt1234
"""

# ---- imports

from collections import defaultdict
from copy import deepcopy
import itertools as it
import pickle
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import StratifiedKFold


from helpers import (get_combo_atn_model,
                     results_boxplot,
                     svm_best_param_lookup,
                     test_atn_linear_model)
from model_classes import ATN_PREDICTORS_DICT, MultivariateSVR

# ---- Variables
OUTER_SPLITS = 5
INNER_SPLITS = 5
REPEATS = 20
OUTER_SEED = 1000
INNER_SEED = 1100

TARGET = 'PACC'
COVARIATES = ['Age', 'Sex', 'HasE4']
STRATIFY_COLUMN = 'CDRBinned'

# ---- Paths
PATH_ADNI_DATA = '../../data/derivatives/adni_harmonized_augmented.csv'
PATH_OASIS_DATA = '../../data/derivatives/oasis_harmonized_augmented.csv'

# ---- Output

OUTPUT = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]
if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

# ---- Read

adni = pd.read_csv(PATH_ADNI_DATA)
adni['Sex'] = np.where(adni['Sex'] == 'Male', 1., 0.)
adni['HasE4'] = adni['HasE4'].astype(float)

oasis = pd.read_csv(PATH_OASIS_DATA)
oasis['Sex'] = np.where(oasis['Sex'] == 'Male', 1., 0.)
oasis['HasE4'] = oasis['HasE4'].astype(float)

# ---- Setup SVM models

amy_columns = list(adni.columns[adni.columns.str.startswith('AV45')])
tau_columns = list(adni.columns[adni.columns.str.startswith('FTP')])
gm_columns = list(adni.columns[adni.columns.str.endswith('VOLUME') & ~ adni.columns.str.contains('BRAAK|META')])
roi_columns = amy_columns + tau_columns + gm_columns

amy_columns += COVARIATES
tau_columns += COVARIATES
gm_columns += COVARIATES
roi_columns += COVARIATES

SVM_MODELS = {'SVM (amyloid)': amy_columns,
              'SVM (tau)': tau_columns,
              'SVM (GM)': gm_columns,
              'SVM (combined)': roi_columns}

SVM_PARAMS = {
    'C': list(2. ** np.arange(-5., 15., 2)),
    'gamma': list(2. ** np.arange(-5., 15., 2)),
    'kernel': ['rbf']}

param_combos = list(it.product(*SVM_PARAMS.values()))
SVM_SEARCH = [dict(zip(SVM_PARAMS.keys(), v)) for v in param_combos]

# ---- Main
results_adni = []
results_oasis = []
save_models = defaultdict(list)

def testing_filter(dataset):
    return dataset.loc[dataset['CDRBinned'] != '0.0', ].copy()

oasis = testing_filter(oasis)

# repeats of nested CV
for r in range(REPEATS):

    outer_cv = StratifiedKFold(n_splits=OUTER_SPLITS, random_state=OUTER_SEED + r, shuffle=True)
    inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, random_state=INNER_SEED + r, shuffle=True)

    # outer CV loop
    for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(adni, adni[STRATIFY_COLUMN])):

        msg = f"REPEAT: {r}, OUTER FOLD: {i}"
        print()
        print(msg)
        print('-' * len(msg))

        outer_train = adni.iloc[outer_train_index, :]
        outer_test = testing_filter(adni.iloc[outer_test_index, :])

        inner_cv_lm_results = []
        inner_cv_svm_results = []

        # inner CV loop
        for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[STRATIFY_COLUMN])):

            print(f'*INNER TRAINING FOLD {j}*')

            inner_train = outer_train.iloc[inner_train_index, :]
            inner_test = testing_filter(outer_train.iloc[inner_test_index, :])

            # testing many ATN models
            for atn, measure_dict in ATN_PREDICTORS_DICT.items():
                for measure_type, model_dict in measure_dict.items():
                    for name, model in model_dict.items():

                        # print(f' - {model}')
                        metrics = test_atn_linear_model(models=model,
                                                        covariates=COVARIATES,
                                                        target=TARGET,
                                                        train_data=inner_train,
                                                        test_data=inner_test)

                        row = {'atn': atn,
                               'measure_type': measure_type,
                               'name': name,
                               'fold': j,
                               **metrics}
                        inner_cv_lm_results.append(row)

            # testing many SVM models
            for svm_name, svm_features in SVM_MODELS.items():
                for c, params in enumerate(SVM_SEARCH):
                    # print(f' - {svm_name} ({params})')
                    model = MultivariateSVR(svm_features, TARGET, **params)
                    model.fit(inner_train)
                    preds = model.predict(inner_test)
                    row = {'name': svm_name,
                           **params,
                           'fold': j,
                           'rmse': mean_squared_error(inner_test[TARGET], preds, squared=False),
                           'r2': r2_score(inner_test[TARGET], preds)}
                    inner_cv_svm_results.append(row)

        # select best ATN model
        inner_cv_lm_results = pd.DataFrame(inner_cv_lm_results)
        lm_model_averages = inner_cv_lm_results.groupby(['atn', 'measure_type', 'name'])['rmse'].agg(mean=np.mean, std=np.std).reset_index()
        best_by_measure = lm_model_averages.groupby(['atn', 'measure_type'])['mean'].idxmin()
        lm_selected_models = lm_model_averages.iloc[best_by_measure]

        # select best SVM model
        inner_cv_svm_results = pd.DataFrame(inner_cv_svm_results)
        svm_model_averages = inner_cv_svm_results.groupby(['name'] + list(SVM_PARAMS.keys()))['rmse'].agg(mean=np.mean, std=np.std).reset_index()
        best_by_params = svm_model_averages.groupby('name')['mean'].idxmin()
        svm_selected_models = svm_model_averages.iloc[best_by_params]

        # develop combinations
        FINAL_ATN_MODELS = {
            'All binary': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'binary', 'binary', 'binary'),
            'Categorical A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'categorical', 'binary', 'binary'),
            'Categorical T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'binary', 'categorical', 'binary'),
            'Categorical N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'binary', 'binary', 'categorical'),
            'All categorical': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'categorical', 'categorical', 'categorical'),
            'Continuous A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'continuous', 'binary', 'binary'),
            'Continuous T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'binary', 'continuous', 'binary'),
            'Continuous N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'binary', 'binary', 'continuous'),
            'All continuous': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS_DICT, 'continuous', 'continuous', 'continuous'),
            }

        print()
        print('*SELECTED MODELS*')
        print(lm_selected_models)
        print(svm_selected_models)

        print()
        print('*OUTER TRAINING*')

        for name, model in FINAL_ATN_MODELS.items():
            print(f' - {name} ({model})')

            # testing on ADNI
            metrics = test_atn_linear_model(models=model,
                                            covariates=COVARIATES,
                                            target=TARGET,
                                            train_data=outer_train,
                                            test_data=outer_test)
            row = {'model': name,
                   'fold': i,
                   'repeat': r,
                   **metrics}
            results_adni.append(row)

            # testing on OASIS
            metrics = test_atn_linear_model(models=model,
                                            covariates=COVARIATES,
                                            target=TARGET,
                                            train_data=outer_train,
                                            test_data=oasis)
            row = {'model': name,
                   'fold': i,
                   'repeat': r,
                   **metrics}
            results_oasis.append(row)


            # save model
            save_models[name].append(deepcopy(model))

        for svm_name, svm_features in SVM_MODELS.items():
            best_params = svm_best_param_lookup(svm_selected_models, svm_name, list(SVM_PARAMS.keys()))
            model = MultivariateSVR(svm_features, TARGET, **best_params)
            model.fit(outer_train)

            print(f' - {svm_name} ({best_params})')

            # test on ADNI
            preds = model.predict(outer_test)
            row = {'model': svm_name,
                   'fold': i,
                   'repeat': r,
                   'rmse': mean_squared_error(outer_test[TARGET], preds, squared=False),
                   'r2': r2_score(outer_test[TARGET], preds)}
            results_adni.append(row)

            # test on OASIS
            preds = model.predict(oasis)
            row = {'model': svm_name,
                   'fold': i,
                   'repeat': r,
                   'rmse': mean_squared_error(oasis[TARGET], preds, squared=False),
                   'r2': r2_score(oasis[TARGET], preds)}
            results_oasis.append(row)

            # save model
            save_models[svm_name].append(deepcopy(model))

results_adni = pd.DataFrame(results_adni)
results_adni.to_csv(os.path.join(OUTPUT, 'adni_results.csv'), index=False)

results_oasis = pd.DataFrame(results_oasis)
results_oasis.to_csv(os.path.join(OUTPUT, 'oasis_results.csv'), index=False)

with open(os.path.join(OUTPUT, 'models.pickle'), 'wb') as f:
    pickle.dump(save_models, f)

#%%

# ---- save plots
plt.rcParams.update({'font.family': 'Arial',
                      'font.size': 15})
n_test = int((1/OUTER_SPLITS) * len(adni))
n_train = len(adni) - n_test

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

# colors
palette = (['gray'] +
           ['#E59C9C'] * 3 +
           ['#ba635d'] +
           ['#AAC3E9'] * 3 +
           ['#7DA1D8'] +
           ['#B3E2AD'] * 3 +
           ['#99C494'])

name = 'adni_rmse_boxplot.png'
title = 'Accuracy (ADNI)'
adni_plot, adni_stats = results_boxplot(results_adni,
                                        save=os.path.join(OUTPUT, name),
                                        title=title,
                                        n_test=n_test,
                                        n_train=n_train,
                                        baseline='All binary',
                                        palette=palette)

name = 'oasis_rmse_boxplot.png'
title = 'Accuracy (OASIS)'
oasis_plot, oasis_stats = results_boxplot(results_oasis,
                                          save=os.path.join(OUTPUT, name),
                                          title=title,
                                          n_test=n_test,
                                          n_train=n_train,
                                          baseline='All binary',
                                          palette=palette)