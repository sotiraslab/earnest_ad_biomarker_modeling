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

import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import mean_squared_error, r2_score

from atn_modeling.atn_predictor_classes import MultivariateSVR
from atn_modeling.atn_predictor_instances import ATN_PREDICTORS, ATN_PREDICTORS_PLUS_CSF
from atn_modeling.helpers import test_atn_linear_model, svm_best_param_lookup, get_combo_atn_model

def experiment_test_all_atn_predictors(dataset, target,
                                       covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                       stratify='CDRBinned',
                                       splits=10,
                                       repeats=10,
                                       seed=0,
                                       savepath=None,
                                       savemodels=None,
                                       with_csf=False):
    
    # hold results
    results = []
    models = defaultdict(list)
    
    start_time = time.time()
    
    predictor_dict = ATN_PREDICTORS_PLUS_CSF if with_csf else ATN_PREDICTORS
    
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
            metrics, _  = test_atn_linear_model(models=[],
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
            for biomarker, variable_dict in predictor_dict.items():
                for variable_type, model_list in variable_dict.items():
                    for model in model_list:
                        metrics, _ = test_atn_linear_model(models=model,
                                                           covariates=covariates,
                                                           target=target,
                                                           train_data=outer_train,
                                                           test_data=outer_test)
                        row = {'biomarker': biomarker,
                               'variable_type': variable_type,
                               'name': model.nickname,
                               'fold': i,
                               'repeat': r,
                               'ntrain': len(outer_train),
                               'ntest': len(outer_test),
                               **metrics}
                        results.append(row)
                        models[model.nickname].append(model)
                        
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
    
    if savemodels:
        print('')
        print(f'Saving models to "{savemodels}"...')
        with open(savemodels, 'wb') as f:
            pickle.dump(models, f)
        print('Done!')
    
    return results_df, models

def experiment_svm(dataset, target,
                   covariates=['Age', 'SexBinary', 'HasE4Binary'],
                   stratify='CDRBinned',
                   search_C_linear=None,
                   search_C_rbf=None,
                   search_kernel=None,
                   search_gamma=None,
                   repeats=10,
                   outer_splits=10,
                   inner_splits=5,
                   outer_seed=0,
                   inner_seed=100,
                   savepath=None,
                   savemodels=None,
                   testing_filter=None):
    
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
    
    if search_C_linear is None:
        search_C_linear = [2 ** -7]

    if search_C_rbf is None:
        search_C_rbf = [2 ** -7]
    
    if search_gamma is None:
        search_gamma = [2 ** -7]
        
    searched_svm_paramters = ['C', 'gamma', 'kernel']
        
    SVM_SEARCH = []
    if 'linear' in search_kernel:
        linear_params = {'C': search_C_linear, 'gamma': [None], 'kernel': ['linear']}
        linear_combos =  list(it.product(*linear_params.values()))
        linear_search = [dict(zip(linear_params.keys(), v)) for v in linear_combos]
        SVM_SEARCH += linear_search
    if 'rbf' in search_kernel:
        rbf_params = {'C': search_C_rbf, 'gamma': search_gamma, 'kernel': ['rbf']}
        rbf_combos =  list(it.product(*rbf_params.values()))
        rbf_search = [dict(zip(rbf_params.keys(), v)) for v in rbf_combos]
        SVM_SEARCH += rbf_search
    
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
            
            if testing_filter is not None:
                outer_test = testing_filter(outer_test)
            
            # inner CV loop
            inner_cv_svm_results = []
            for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[stratify])):
    
                print(f'[{str(dt.datetime.now())}] *INNER TRAINING FOLD {j}*')
    
                inner_train = outer_train.iloc[inner_train_index, :]
                inner_test = outer_train.iloc[inner_test_index, :]
                
                if testing_filter is not None:
                    inner_test = testing_filter(inner_test)
                
                # testing many SVM models
                
                for svm_name, svm_features in SVM_MODELS.items():
                    for c, params in enumerate(SVM_SEARCH):
                        # uncomment to print every model trained
                        # much more verbose!
                        # print(f' - {svm_name} ({params})')
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
            svm_model_averages = inner_cv_svm_results.groupby(['name'] + searched_svm_paramters)['rmse'].agg(mean="mean", std="std").reset_index()
            best_by_params = svm_model_averages.groupby('name')['mean'].idxmin()
            svm_selected_models = svm_model_averages.iloc[best_by_params]
            
            print()
            print(' SELECTED MODELS*')
            print(svm_selected_models)
    
            print()
            print(f'[{str(dt.datetime.now())}] *OUTER TRAINING*')
            for svm_name, svm_features in SVM_MODELS.items():
                best_params = svm_best_param_lookup(svm_selected_models, svm_name, searched_svm_paramters)
                model = MultivariateSVR(svm_features, target, **best_params)
                model.fit(outer_train)
    
                print(f' - {svm_name} ({best_params}) [{str(dt.datetime.now())}]')
    
                # test on ADNI
                preds = model.predict(outer_test)
                row = {'model': svm_name,
                       **best_params, 
                       'biomarker': BIOMARKERS[svm_name],
                       'fold': i,
                       'repeat': r,
                       'ntrain': len(outer_train),
                       'ntest': len(outer_test),
                       'rmse': mean_squared_error(outer_test[target], preds, squared=False),
                       'r2': r2_score(outer_test[target], preds)}
                results.append(row)
                print(f'   RMSE={row["rmse"]}')

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

def experiment_combo_atn_vs_baseline(dataset, target,
                                     covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                     stratify='CDRBinned',
                                     repeats=10,
                                     outer_splits=10,
                                     inner_splits=5,
                                     outer_seed=0,
                                     inner_seed=100,
                                     savepath=None,
                                     savemodels=None,
                                     savelms=None,
                                     testing_filter=None):
    
    results = []
    models = defaultdict(list)
    lms = defaultdict(list)
    
    # repeats of nested CV
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
            
            if testing_filter is not None:
                outer_test = testing_filter(outer_test)
    
            inner_cv_lm_results = []
    
            # inner CV loop
            for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[stratify])):
    
                print(f'[{str(dt.datetime.now())}] *INNER TRAINING FOLD {j}*')
    
                inner_train = outer_train.iloc[inner_train_index, :]
                inner_test = outer_train.iloc[inner_test_index, :]
                
                if testing_filter is not None:
                    inner_test = testing_filter(inner_test)
    
                # testing many ATN models
                for biomarker, variable_dict in ATN_PREDICTORS.items():
                    for variable_type, model_list in variable_dict.items():
                        for model in model_list:
                            metrics, _ = test_atn_linear_model(models=model,
                                                               covariates=covariates,
                                                               target=target,
                                                               train_data=inner_train,
                                                               test_data=inner_test)
                            row = {'biomarker': biomarker,
                                   'variable_type': variable_type,
                                   'name': model.nickname,
                                   'fold': i,
                                   'repeat': r,
                                   **metrics}
                            inner_cv_lm_results.append(row)
    
            # select best ATN model
            inner_cv_lm_results = pd.DataFrame(inner_cv_lm_results)
            lm_model_averages = inner_cv_lm_results.groupby(['biomarker', 'variable_type', 'name'])['rmse'].agg(mean='mean', std='std').reset_index()
            best_by_measure = lm_model_averages.groupby(['biomarker', 'variable_type'])['mean'].idxmin()
            lm_selected_models = lm_model_averages.iloc[best_by_measure]
    
            # develop combinations
            FINAL_ATN_MODELS = {
                'Baseline': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, None),
                'Binary A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', None, None),
                'Binary T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'binary', None),
                'Binary T [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'binary', None, taupvc=True),
                'Binary N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'binary'),
                'All binary': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'binary'),
                'All binary [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'binary', taupvc=True),
                'Categorical A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', None, None),
                'Categorical T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'categorical', None),
                'Categorical T [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'categorical', None, taupvc=True),
                'Categorical N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'categorical'),
                'All categorical': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'categorical', 'categorical'),
                'All categorical [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'categorical', 'categorical', taupvc=True),
                'Continuous A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', None, None),
                'Continuous T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'continuous', None),
                'Continuous T [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, 'continuous', None, taupvc=True),
                'Continuous N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, None, None, 'continuous'),
                'All continuous': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'continuous', 'continuous'),
                'All continuous [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'continuous', 'continuous', taupvc=True),
                }
    
            print()
            print('[{str(dt.datetime.now())}] *OUTER TRAINING*')
    
            for name, model in FINAL_ATN_MODELS.items():
                print(f' - {name} ({[m.nickname for m in model]})')
    
                # testing on ADNI
                metrics, lm = test_atn_linear_model(models=model,
                                                    covariates=covariates,
                                                    target=target,
                                                    train_data=outer_train,
                                                    test_data=outer_test)
                row = {'model': name,
                       'fold': i,
                       'repeat': r,
                       'ntrain': len(outer_train),
                       'ntest': len(outer_test),
                       **metrics}
                results.append(row)
    
                # save model
                models[name].append(deepcopy(model))
                lms[name].append(deepcopy(lm))

    
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
        
    if savelms:
        print('')
        print(f'Saving linear models to "{savelms}"...')
        with open(savelms, 'wb') as f:
            pickle.dump(lms, f)
        print('Done!')
            
    return results_df, models, lms

def experiment_combo_atn_vs_binary(dataset, target,
                                   covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                   stratify='CDRBinned',
                                   repeats=10,
                                   outer_splits=10,
                                   inner_splits=5,
                                   outer_seed=0,
                                   inner_seed=100,
                                   savepath=None,
                                   savemodels=None,
                                   savelms=None,
                                   testing_filter=None):
    
    results = []
    models = defaultdict(list)
    lms = defaultdict(list)
    
    # repeats of nested CV
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
            
            if testing_filter is not None:
                outer_test = testing_filter(outer_test)
    
            inner_cv_lm_results = []
    
            # inner CV loop
            for j, (inner_train_index, inner_test_index) in enumerate(inner_cv.split(outer_train, outer_train[stratify])):
    
                print(f'[{str(dt.datetime.now())}] *INNER TRAINING FOLD {j}*')
    
                inner_train = outer_train.iloc[inner_train_index, :]
                inner_test = outer_train.iloc[inner_test_index, :]
                
                if testing_filter is not None:
                    inner_test = testing_filter(inner_test)
    
                # testing many ATN models
                for biomarker, variable_dict in ATN_PREDICTORS.items():
                    for variable_type, model_list in variable_dict.items():
                        for model in model_list:
                            metrics, _ = test_atn_linear_model(models=model,
                                                               covariates=covariates,
                                                               target=target,
                                                               train_data=inner_train,
                                                               test_data=inner_test)
                            row = {'biomarker': biomarker,
                                   'variable_type': variable_type,
                                   'name': model.nickname,
                                   'fold': i,
                                   'repeat': r,
                                   **metrics}
                            inner_cv_lm_results.append(row)
    
            # select best ATN model
            inner_cv_lm_results = pd.DataFrame(inner_cv_lm_results)
            lm_model_averages = inner_cv_lm_results.groupby(['biomarker', 'variable_type', 'name'])['rmse'].agg(mean='mean', std='std').reset_index()
            best_by_measure = lm_model_averages.groupby(['biomarker', 'variable_type'])['mean'].idxmin()
            lm_selected_models = lm_model_averages.iloc[best_by_measure]
    
            # develop combinations
            FINAL_ATN_MODELS = {
            'All binary': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'binary'),
            'Categorical A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'binary', 'binary'),
            'Categorical T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'categorical', 'binary'),
            'Categorical N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'categorical'),
            'All categorical': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'categorical', 'categorical'),
            'Continuous A': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'binary', 'binary'),
            'Continuous T': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'continuous', 'binary'),
            'Continuous N': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'continuous'),
            'All continuous': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'continuous', 'continuous'),
            'All binary [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'binary', taupvc=True),
            'Categorical A [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'binary', 'binary', taupvc=True),
            'Categorical T [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'categorical', 'binary', taupvc=True),
            'Categorical N [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'categorical', taupvc=True),
            'All categorical [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'categorical', 'categorical', 'categorical', taupvc=True),
            'Continuous A [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'binary', 'binary', taupvc=True),
            'Continuous T [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'continuous', 'binary', taupvc=True),
            'Continuous N [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'binary', 'binary', 'continuous', taupvc=True),
            'All continuous [PVC]': get_combo_atn_model(lm_selected_models, ATN_PREDICTORS, 'continuous', 'continuous', 'continuous', taupvc=True),
            }
    
            print()
            print('[{str(dt.datetime.now())}] *OUTER TRAINING*')
    
            for name, model in FINAL_ATN_MODELS.items():
                print(f' - {name} ({[m.nickname for m in model]})')
    
                # testing on ADNI
                metrics, lm = test_atn_linear_model(models=model,
                                                    covariates=covariates,
                                                    target=target,
                                                    train_data=outer_train,
                                                    test_data=outer_test)
                row = {'model': name,
                       'fold': i,
                       'repeat': r,
                       'ntrain': len(outer_train),
                       'ntest': len(outer_test),
                       **metrics}
                results.append(row)
    
                # save model
                models[name].append(deepcopy(model))
                lms[name].append(deepcopy(lm))

    
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
        
    if savelms:
        print('')
        print(f'Saving linear models to "{savelms}"...')
        with open(savelms, 'wb') as f:
            pickle.dump(lms, f)
        print('Done!')
            
    return results_df, models, lms