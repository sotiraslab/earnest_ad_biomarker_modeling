#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:16:00 2024

@author: tom.earnest
"""

import time

import pandas as pd
from sklearn.model_selection import StratifiedKFold

from atn_predictor_instances import ATN_PREDICTORS
from helpers import test_atn_linear_model

def experiment_test_all_atn_predictors(dataset, target,
                                       covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                       stratify='CDRBinned',
                                       splits=5,
                                       repeats=5,
                                       seed=0):
    
    # hold results
    results = []
    
    # repeats of cross validation
    start_time = time.time()
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
                        
            # test baseline model (no ATN info)
            metrics = test_atn_linear_model(models=[],
                                            covariates=covariates,
                                            target=target,
                                            train_data=outer_train,
                                            test_data=outer_test)
            row = {'biomarker': 'NA',
                   'variable_type': 'NA',
                   'name': 'baseline',
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
    return results_df
