#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import argparse
import os

import numpy as np
import pandas as pd

from experiments import experiment_svm
from helpers import results_boxplot

def parse():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--short', action='store_true')
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--noplot', action='store_true')
    
    return parser.parse_args()

def main(rerun=False, replot=True, short=False):
    
    # setup output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    if short:
        experiment_name += '_short'
    this_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    # check existing output
    results_path = os.path.join(output_folder, 'results.csv')
    models_path = os.path.join(output_folder, 'models.pickle')
    results_exist = os.path.exists(results_path)
        
    # EXPERIMENT PARAMETERS
    dataset = pd.read_csv('../../outputs/maindata/maindata.csv')
    target = 'PHC_GLOBAL'
    covariates=['Age', 'SexBinary', 'HasE4Binary']
    stratify='CDRBinned'
    search_C_linear = [2 ** -7] if short else list(2. ** np.arange(-10, 1, 1))
    search_C_rbf = [2] if short else list(2. ** np.arange(-5, 17, 2))
    search_gamma = [2 ** -11] if short else list(2. ** np.arange(-15, 5, 2))
    search_kernel = ['linear', 'rbf']
    repeats=10
    outer_splits=10
    inner_splits=5
    outer_seed=0
    inner_seed=100
    
    # run
    if rerun or (not results_exist):
        _ = experiment_svm(dataset,
                           target=target,
                           covariates=covariates,
                           stratify=stratify,
                           repeats=repeats, 
                           outer_splits=outer_splits,
                           inner_splits=inner_splits,
                           outer_seed=outer_seed,
                           inner_seed=inner_seed,
                           search_C_linear=search_C_linear,
                           search_C_rbf=search_C_rbf,
                           search_kernel=search_kernel,
                           search_gamma=search_gamma,
                           savepath=results_path,
                           savemodels=models_path)
    results = pd.read_csv(results_path)
        
    if replot:
        
        # general resources
        plot_path = os.path.join(output_folder, 'boxplot.svg')
        colors = {'Amyloid SVM': '#882255',
                  'Tau SVM': '#117733',
                  'Tau SVM [PVC]': '#117733',
                  'GM SVM': '#332288',
                  'ATN SVM': '#DDCC77',
                  'ATN SVM [PVC]': '#DDCC77',
                  None: 'Gray'}
        
        # 1. grouped by biomarker
        data = results.copy()
        order = data['model'].unique()
        group = data.groupby('model').agg({'model': 'first', 'rmse': 'mean'})
        group['color'] = group['model'].map(colors)
        group = group.loc[order, :]
        
        fig, _ = results_boxplot(data, groupby='model', baseline=None,
                                 stats=False, palette=group['color'],
                                 save=plot_path)

    
if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    kwargs['replot'] = not kwargs['noplot']
    del kwargs['noplot']
    main(**kwargs)