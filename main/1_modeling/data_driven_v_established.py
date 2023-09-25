#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:22:26 2023

@author: earnestt1234
"""

# ---- imports

from collections import defaultdict
from copy import deepcopy
import pickle
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold


from helpers import (get_combo_atn_model,
                     results_boxplot,
                     test_atn_linear_model)
from model_classes import BINARY_DATA_DRIVEN, BINARY_ESTABLISHED

# ---- Variables
OUTER_SPLITS = 5
INNER_SPLITS = 5
REPEATS = 1
OUTER_SEED = 1600
INNER_SEED = 1700

TARGET = 'PACC'
COVARIATES = ['Age', 'Sex', 'HasE4']
STRATIFY_COLUMN = 'CDRBinned'

# ---- Paths
PATH_ADNI_DATA = '../../data/derivatives/adni_base_table.csv'

# ---- Output

OUTPUT = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]
if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

# ---- Read

adni = pd.read_csv(PATH_ADNI_DATA)
adni['Sex'] = np.where(adni['Sex'] == 'Male', 1., 0.)
adni['HasE4'] = adni['HasE4'].astype(float)

# ---- Main
results_adni = []
save_models = defaultdict(list)

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
        outer_test = adni.iloc[outer_test_index, :]

        inner_cv_lm_results = []
        inner_cv_svm_results = []

        # inner CV loop
        for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[STRATIFY_COLUMN])):

            print(f'*INNER TRAINING FOLD {j}*')

            inner_train = outer_train.iloc[inner_train_index, :]
            inner_test = outer_train.iloc[inner_test_index, :]

            # testing many ATN models
            for atn, measure_dict in BINARY_DATA_DRIVEN.items():
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
                               'learning_type': 'datadriven',
                               'name': name,
                               'fold': j,
                               **metrics}
                        inner_cv_lm_results.append(row)
                        
            for atn, measure_dict in BINARY_ESTABLISHED.items():
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
                               'learning_type': 'established',
                               'name': name,
                               'fold': j,
                               **metrics}
                        inner_cv_lm_results.append(row)


        # select best ATN model
        inner_cv_lm_results = pd.DataFrame(inner_cv_lm_results)
        lm_model_averages = inner_cv_lm_results.groupby(['atn', 'measure_type', 'learning_type', 'name'])['rmse'].agg(mean=np.mean, std=np.std).reset_index()
        best_by_measure = lm_model_averages.groupby(['atn', 'measure_type', 'learning_type'])['mean'].idxmin()
        lm_selected_models = lm_model_averages.iloc[best_by_measure]
        
        lm_datadriven = lm_selected_models[lm_selected_models['learning_type'] == 'datadriven']
        lm_established = lm_selected_models[lm_selected_models['learning_type'] == 'established']

        # develop combinations
        FINAL_ATN_MODELS = {
            'Baseline': get_combo_atn_model(lm_selected_models, BINARY_DATA_DRIVEN, None, None, None),
            'Literature A': get_combo_atn_model(lm_established, BINARY_ESTABLISHED, 'binary', None, None),
            'Literature T': get_combo_atn_model(lm_established, BINARY_ESTABLISHED, None, 'binary', None),
            'Literature A+T': get_combo_atn_model(lm_established, BINARY_ESTABLISHED, 'binary', 'binary', None),
            'Data-driven A': get_combo_atn_model(lm_datadriven, BINARY_DATA_DRIVEN, 'binary', None, None),
            'Data-driven T': get_combo_atn_model(lm_datadriven, BINARY_DATA_DRIVEN, None, 'binary', None),
            'Data-driven A+T': get_combo_atn_model(lm_datadriven, BINARY_DATA_DRIVEN, 'binary', 'binary', None),
            }

        print()
        print('*SELECTED MODELS*')
        print(lm_selected_models)

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

            # save model
            save_models[name].append(deepcopy(model))

results_adni = pd.DataFrame(results_adni)
results_adni.to_csv(os.path.join(OUTPUT, 'adni_results.csv'), index=False)

with open(os.path.join(OUTPUT, 'models.pickle'), 'wb') as f:
    pickle.dump(save_models, f)

# ---- save plots
plt.rcParams.update({'font.family': 'Arial',
                      'font.size': 15})
n_test = int((1/OUTER_SPLITS) * len(adni))
n_train = len(adni) - n_test

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

# colors
palette = (['gray'] +
            ['#DDCBEF'] * 2 +
            ['#AA6CE9'] +
            ['#C1ECDB'] * 2 +
            ['#55C79A'])

name = 'adni_rmse_boxplot.png'
title = 'Accuracy (ADNI)'
adni_plot, adni_stats = results_boxplot(results_adni,
                                        save=os.path.join(OUTPUT, name),
                                        title=title,
                                        n_test=n_test,
                                        n_train=n_train,
                                        baseline='Baseline',
                                        palette=palette)
