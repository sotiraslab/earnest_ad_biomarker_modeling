#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:16:00 2024

@author: tom.earnest
"""

import time

from sklearn.metrics import mean_squared_error, r2_score
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
            print(f'Time: {seconds}s')
                
                        
    results_df = pd.DataFrame(results)
    return results_df

import pandas as pd
dataset = pd.read_csv('../../outputs/datasets/maindata.csv')
results = experiment_test_all_atn_predictors(dataset, target='PHC_GLOBAL',
                                             repeats=10, splits=10)

#%%
import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

transposed = results.pivot(index=['repeat', 'fold'], columns='name', values='rmse')
transposed = transposed.loc[:, ~ transposed.columns.str.contains('PVC')]

fig, ax = plt.subplots(figsize=(6, 12))

for i in range(transposed.shape[1]):
    values = transposed.iloc[:, i]
    kde = gaussian_kde(values)
    x = np.linspace(transposed.min(axis=None), transposed.max(axis=None), 100)
    y = kde(x)
    y /= y.max()
    y += i
    # ax.plot(x, y)
    ax.fill_between(x, y1=y, y2=i, alpha=.3)
    
ax.set_yticks(range(transposed.shape[1]))
ax.set_yticklabels(transposed.columns)



