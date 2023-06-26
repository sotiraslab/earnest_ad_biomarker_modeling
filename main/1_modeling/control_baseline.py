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
from model_classes import (BinaryGMM,
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
OUTER_SEED = 0
INNER_SEED = 100

TARGET = 'PACC.ADNI'
COVARIATES = ['Age', 'Sex', 'HasE4']
STRATIFY_COLUMN = 'CDRBinned'

# ---- Paths
PATH_DATA = '../../data/derivatives/adni_base_table.csv'


# ---- Output

OUTPUT = os.path.splitext(os.path.basename(os.path.abspath(__file__)))[0]
if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)

# ---- Read

df = pd.read_csv(PATH_DATA)
df['Sex'] = np.where(df['Sex'] == 'Male', 1., 0.)
df['HasE4'] = df['HasE4'].astype(float)

# ---- ATN predeictors

ATN_PREDICTORS = {
    'amyloid': {
        'binary': {
            'composite_wm_1.11': BinaryManual('SUMMARYSUVR_WHOLECEREBNORM', 1.11),
            'centiloid_20': BinaryManual('Centiloid', 20)},
        'categorical': {
            'composite_wm_quantiles': Quantiles('SUMMARYSUVR_WHOLECEREBNORM'),
            'centiloid_quantiles': Quantiles('Centiloid')},
        'continuous': {
            'composite_wm': Continuous('SUMMARYSUVR_WHOLECEREBNORM'),
            'centiloid': Continuous('Centiloid')}},
    'tau': {
        'binary': {
            'mtt_gmm': BinaryGMM('META_TEMPORAL_SUVR'),
            'mtt_z2.0': BinaryZScore('META_TEMPORAL_SUVR', 'Control', zcutoff=2.0),
            'mtt_z2.5': BinaryZScore('META_TEMPORAL_SUVR', 'Control', zcutoff=2.5),
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
            'hipp_z2.0': BinaryZScore('HIPPOCAMPUS_VOLUME', control_col='Control', zcutoff=-2.0, greater=False),
            'hipp_z2.5': BinaryZScore('HIPPOCAMPUS_VOLUME', control_col='Control', zcutoff=-2.5, greater=False),
            'mttvol_z2.0': BinaryZScore('META_TEMPORAL_VOLUME', control_col='Control', zcutoff=-2.0, greater=False),
            'mttvol_z2.5': BinaryZScore('META_TEMPORAL_VOLUME', control_col='Control', zcutoff=-2.5, greater=False)},
        'categorical': {
            'hipp_quantiles': Quantiles('HIPPOCAMPUS_VOLUME'),
            'mttvol_quantiles': Quantiles('META_TEMPORAL_VOLUME')},
        'continuous': {
            'hipp': Continuous('HIPPOCAMPUS_VOLUME'),
            'mttvol': Continuous('META_TEMPORAL_VOLUME')}},
    }

# ---- Setup models

amy_columns = list(df.columns[df.columns.str.startswith('AV45')])
tau_columns = list(df.columns[df.columns.str.startswith('FTP')])
gm_columns = list(df.columns[df.columns.str.endswith('VOLUME') & ~ df.columns.str.contains('BRAAK|META')])
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
    'C': [0.01, 0.1, 1, 10, 100],
    'epsilon': [0.01, .1, .5, 1, 2],
    'kernel': ['rbf', 'linear']}

param_combos = list(it.product(*SVM_PARAMS.values()))
SVM_SEARCH = [dict(zip(SVM_PARAMS.keys(), v)) for v in param_combos]

# ---- Main
results = []
save_models = defaultdict(list)

# repeats of nested CV
for r in range(REPEATS):

    outer_cv = StratifiedKFold(n_splits=OUTER_SPLITS, random_state=OUTER_SEED + r, shuffle=True)
    inner_cv = StratifiedKFold(n_splits=INNER_SPLITS, random_state=INNER_SEED + r, shuffle=True)

    # outer CV loop
    for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(df, df[STRATIFY_COLUMN])):
        
        msg = f"REPEAT: {r}, OUTER FOLD: {i}"
        print()
        print(msg)
        print('-' * len(msg))
        
        outer_train = df.iloc[outer_train_index, :]
        outer_test = df.iloc[outer_test_index, :]

        inner_cv_lm_results = []
        inner_cv_svm_results = []

        # inner CV loop
        for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[STRATIFY_COLUMN])):
            
            print(f'*INNER TRAINING FOLD {j}*')
            
            inner_train = outer_train.iloc[inner_train_index, :]
            inner_test = outer_train.iloc[inner_test_index, :]

            # testing many ATN models
            for atn, measure_dict in ATN_PREDICTORS.items():
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
            'Baseline': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, None),
            'Binary A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', None, None),
            'Binary T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'binary', None),
            'Binary N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'binary'),
            'All binary': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'binary'),
            'Categorical A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', None, None),
            'Categorical T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'categorical', None),
            'Categorical N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'categorical'),
            'All categorical': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'categorical', 'categorical'),
            'Continuous A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', None, None),
            'Continuous T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'continuous', None),
            'Continuous N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'continuous'),
            'All continuous': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'continuous', 'continuous'),
            }
        
        print()
        print('*SELECTED MODELS*')
        print(lm_selected_models)
        print(svm_selected_models)
        
        print()
        print('*OUTER TRAINING*')

        for name, model in FINAL_ATN_MODELS.items():
            print(f' - {name} ({model})')
            metrics = test_atn_linear_model(models=model,
                                            covariates=COVARIATES,
                                            target=TARGET, 
                                            train_data=outer_train, 
                                            test_data=outer_test)
            row = {'model': name,
                   'fold': i,
                   'repeat': r,
                   **metrics}
            results.append(row)
            save_models[name].append(deepcopy(model))

        for svm_name, svm_features in SVM_MODELS.items():
            best_params = svm_best_param_lookup(svm_selected_models, svm_name, list(SVM_PARAMS.keys()))
            print(f' - {svm_name} ({best_params})')
            model = MultivariateSVR(svm_features, TARGET, **best_params)
            model.fit(outer_train)
            preds = model.predict(outer_test)
            row = {'model': name,
                   'fold': i,
                   'repeat': r,
                   'rmse': mean_squared_error(outer_test[TARGET], preds, squared=False),
                   'r2': r2_score(outer_test[TARGET], preds)}
            results.append(row)
            save_models[svm_name].append(deepcopy(model))

    break

results = pd.DataFrame(results)
results.to_csv(os.path.join(OUTPUT, 'outer_cv_results.csv'), index=False)

with open(os.path.join(OUTPUT, 'models.pickle'), 'wb') as f:
    pickle.dump(save_models, f)

# ---- save plot
plt.rcParams.update({'font.family': 'Arial',
                      'font.size': 15})
n_test = int((1/OUTER_SPLITS) * len(df))
n_train = len(df) - n_test

if not os.path.isdir(OUTPUT):
    os.mkdir(OUTPUT)
    
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

name = 'rmse_boxplot.png'
title = 'Accuracy (entire test set)'
results_boxplot(results,
                save=os.path.join(OUTPUT, name),
                title=title,
                n_test=n_test,
                n_train=n_train,
                baseline='Baseline',
                palette=palette)