#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import argparse
import os

import pandas as pd

from experiments import experiment_combo_atn_vs_binary
from helpers import results_boxplot

def parse():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--noplot', action='store_true')
    
    return parser.parse_args()

def main(rerun=False, plot=True):
    
    # setup output
    experiment_name = os.path.splitext(os.path.basename(__file__))[0]
    this_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', experiment_name))
    
    if not os.path.isdir(output_folder):
        os.mkdir(output_folder)
        
    # check existing output
    results_path = os.path.join(output_folder, 'results.csv')
    models_path = os.path.join(output_folder, 'models.pickle')
    lms_path = os.path.join(output_folder, 'lm.pickle')
    results_exist = os.path.exists(results_path)
        
    # run experiment
    if rerun or (not results_exist):
        
        # EXPERIMENT PARAMETERS
        
        # NOTE
        # seeds here should be same as exp1
        # for combining of results
        
        dataset = pd.read_csv('../../outputs/maindata/maindata.csv')
        target = 'PHC_EXF'
        covariates=['Age', 'SexBinary', 'HasE4Binary']
        stratify='CDRBinned'
        repeats=10
        outer_splits=10
        inner_splits=5
        outer_seed=0
        inner_seed=100
        
        dataset = pd.read_csv('../../outputs/maindata/maindata.csv')
        _ = experiment_combo_atn_vs_binary(dataset,
                                           target=target,
                                           covariates=covariates,
                                           stratify=stratify,
                                           repeats=repeats, 
                                           outer_splits=outer_splits,
                                           inner_splits=inner_splits,
                                           outer_seed=outer_seed,
                                           inner_seed=inner_seed,
                                           savepath=results_path,
                                           savemodels=models_path,
                                           savelms=lms_path)
        
    results = pd.read_csv(results_path)
        
    if plot:
        
        # general resources
        plot_path = os.path.join(output_folder, 'boxplot.svg')
        
        # plot
        data = results.loc[~results['model'].str.contains('PVC'), :]
        
        palette = (['gray'] +
            ['#E59C9C'] * 3 +
            ['#ba635d'] +
            ['#AAC3E9'] * 3 +
            ['#7DA1D8'] +
            ['#B3E2AD'] * 3 +
            ['#99C494'])
        
        fig, _ = results_boxplot(data, groupby='model', baseline='All binary',
                                 stats=False, palette=palette,
                                 save=plot_path)
    
if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    kwargs['plot'] = not kwargs['noplot']
    del kwargs['noplot']
    main(**kwargs)