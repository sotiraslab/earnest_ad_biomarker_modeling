#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:16:00 2024

@author: tom.earnest
"""

from collections import defaultdict
from copy import deepcopy
import datetime as dt
import itertools as it
import pickle
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import mean_squared_error, r2_score

from atn_predictor_classes import MultivariateSVR
from atn_predictor_instances import ATN_PREDICTORS
from helpers import test_atn_linear_model, svm_best_param_lookup

def experiment_test_all_atn_predictors(dataset, target,
                                       covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                       stratify='CDRBinned',
                                       splits=10,
                                       repeats=10,
                                       seed=0,
                                       savepath=None):
    
    # hold results
    results = []
    
    start_time = time.time()
    
    # repeats of cross validation
    for r in range(repeats):
        
        cv = StratifiedKFold(n_splits=splits, random_state=r+seed, shuffle=True)
        for i, (outer_train_index, outer_test_index) in enumerate(cv.split(dataset, dataset[stratify])):
            
            msg = f"REPEAT: {r}, FOLD: {i}"
            print()
            print(msg)
            print('-' * len(msg))
            start = time.time()
            
            outer_train = dataset.iloc[outer_train_index, :]
            outer_test = dataset.iloc[outer_test_index, :]
            
            # test baseline model (no ATN info)
            metrics = test_atn_linear_model(models=[],
                                            covariates=covariates,
                                            target=target,
                                            train_data=outer_train,
                                            test_data=outer_test)
            row = {'biomarker': 'NA',
                   'variable_type': 'NA',
                   'name': 'Baseline',
                   'fold': i,
                   'repeat': r,
                   **metrics}
            results.append(row)
            
            # test ATN models
            for biomarker, variable_dict in ATN_PREDICTORS.items():
                for variable_type, model_list in variable_dict.items():
                    for model in model_list:
                        metrics = test_atn_linear_model(models=model,
                                                        covariates=covariates,
                                                        target=target,
                                                        train_data=outer_train,
                                                        test_data=outer_test)
                        row = {'biomarker': biomarker,
                               'variable_type': variable_type,
                               'name': model.nickname,
                               'fold': i,
                               'repeat': r,
                               **metrics}
                        results.append(row)   
                        
            end = time.time()
            seconds = round(end - start, 2)
            elapsed = round(end - start_time, 2)
            print(f'Time: {seconds}s')
            print(f'Elapsed: {elapsed}s')
                        
    results_df = pd.DataFrame(results)
    
    if savepath:
        print('')
        print(f'Saving results to "{savepath}"...')
        results_df.to_csv(savepath, index=False)
        print('Done!')
    
    return results_df

def experiment_svm(dataset, target,
                   covariates=['Age', 'SexBinary', 'HasE4Binary'],
                   stratify='CDRBinned',
                   search_C=None,
                   search_kernel=None,
                   repeats=10,
                   outer_splits=10,
                   inner_splits=5,
                   outer_seed=0,
                   inner_seed=100,
                   savepath=None,
                   savemodels=None):
    
    # setup columns for SVM
    amy_columns = list(dataset.columns[dataset.columns.str.startswith('AV45') & ~dataset.columns.str.contains('TOT')])
    tau_columns = list(dataset.columns[dataset.columns.str.startswith('FTP_') & ~dataset.columns.str.contains('TOT')])
    taupvc_columns = list(dataset.columns[dataset.columns.str.startswith('FTPPVC') & ~dataset.columns.str.contains('TOT')])
    gm_columns = list(dataset.columns[dataset.columns.str.endswith('VOLUME') & ~ dataset.columns.str.contains('BRAAK|META')])
    roi_columns = amy_columns + tau_columns + gm_columns
    roipvc_columns = amy_columns + taupvc_columns + gm_columns
    
    amy_columns += covariates
    tau_columns += covariates
    gm_columns += covariates
    roi_columns += covariates
    roipvc_columns += covariates
    
    assert len(set([len(x) for x in [amy_columns, tau_columns, gm_columns]])) == 1, 'Different # of columns for SVM'

    SVM_MODELS = {'Amyloid SVM': amy_columns,
                  'Tau SVM': tau_columns,
                  'Tau SVM [PVC]': taupvc_columns,
                  'GM SVM': gm_columns,
                  'ATN SVM': roi_columns,
                  'ATN SVM [PVC]': roipvc_columns}
    
    BIOMARKERS = {'Amyloid SVM': 'amyloid',
                  'Tau SVM': 'tau',
                  'Tau SVM [PVC]': 'taupvc',
                  'GM SVM': 'neurodegeneration',
                  'ATN SVM': 'multimodal',
                  'ATN SVM [PVC]': 'multimodal'}
    
    
    if search_C is None:
        search_C = list(2. ** np.arange(-9, 1, 1))
    
    if search_kernel is None:
        search_kernel = ['linear']
        
    SVM_PARAMS = {'C': search_C, 'kernel': search_kernel}
    param_combos = list(it.product(*SVM_PARAMS.values()))
    SVM_SEARCH = [dict(zip(SVM_PARAMS.keys(), v)) for v in param_combos]
    
    # this happens a lot with the linear SVM for certain values of C
    # some models are converging, however
    warnings.filterwarnings('ignore', message='Liblinear failed to converge')
    
    # repeats of cross validation
    results = []
    models = defaultdict(list)
    for r in range(repeats):
        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)
        inner_cv = StratifiedKFold(n_splits=inner_splits, random_state=inner_seed + r, shuffle=True)
        
        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(dataset, dataset[stratify])):
            msg = f"[{str(dt.datetime.now())}] REPEAT: {r}, OUTER FOLD: {i}"
            print()
            print(msg)
            print('-' * len(msg))
            
            outer_train = dataset.iloc[outer_train_index, :]
            outer_test = dataset.iloc[outer_test_index, :]       
            
            # inner CV loop
            inner_cv_svm_results = []
            for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[stratify])):
    
                print(f'[{str(dt.datetime.now())}] *INNER TRAINING FOLD {j}*')
    
                inner_train = outer_train.iloc[inner_train_index, :]
                inner_test = outer_train.iloc[inner_test_index, :]
                
                # testing many SVM models
                for svm_name, svm_features in SVM_MODELS.items():
                    for c, params in enumerate(SVM_SEARCH):
                        print(f' - {svm_name} ({params})')
                        model = MultivariateSVR(svm_features, target, **params)
                        model.fit(inner_train)
                        preds = model.predict(inner_test)
                        row = {'name': svm_name,
                               **params,
                               'repeat': r,
                               'fold': j,
                               'rmse': mean_squared_error(inner_test[target], preds, squared=False),
                               'r2': r2_score(inner_test[target], preds)}
                        inner_cv_svm_results.append(row)
                
            # select best SVM model
            inner_cv_svm_results = pd.DataFrame(inner_cv_svm_results)
            svm_model_averages = inner_cv_svm_results.groupby(['name'] + list(SVM_PARAMS.keys()))['rmse'].agg(mean="mean", std="std").reset_index()
            best_by_params = svm_model_averages.groupby('name')['mean'].idxmin()
            svm_selected_models = svm_model_averages.iloc[best_by_params]
            
            print()
            print(' SELECTED MODELS*')
            print(svm_selected_models)
    
            print()
            print('[{str(dt.datetime.now())}] *OUTER TRAINING*')
            for svm_name, svm_features in SVM_MODELS.items():
                best_params = svm_best_param_lookup(svm_selected_models, svm_name, list(SVM_PARAMS.keys()))
                model = MultivariateSVR(svm_features, target, **best_params)
                model.fit(outer_train)
    
                print(f' - {svm_name} ({best_params}) [{str(dt.datetime.now())}]')
    
                # test on ADNI
                preds = model.predict(outer_test)
                row = {'model': svm_name,
                       'biomarker': BIOMARKERS[svm_name],
                       'fold': i,
                       'repeat': r,
                       'rmse': mean_squared_error(outer_test[target], preds, squared=False),
                       'r2': r2_score(outer_test[target], preds)}
                results.append(row)

                # save model
                models[svm_name].append(deepcopy(model))
            
    results_df = pd.DataFrame(results)
    
    if savepath:
        print('')
        print(f'Saving results to "{savepath}"...')
        results_df.to_csv(savepath, index=False)
        print('Done!')
        
    if savemodels:
        print('')
        print(f'Saving models to "{savemodels}"...')
        with open(savemodels, 'wb') as f:
            pickle.dump(models, f)
        print('Done!')
            
    return results_df, models
        