#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""
import os

import pandas as pd

from atn_modeling.experiments import experiment_test_all_atn_predictors
from common import parse, setup_output

def main(short=False):
    
    output_folder = setup_output(__file__)
        
    # parameters
    dataset=pd.read_csv('datasets/maindata.csv')
    target='PHC_GLOBAL'
    covariates=['Age', 'SexBinary', 'HasE4Binary']
    stratify='CDRBinned'
    repeats=10
    splits=10
    seed=0
    savepath=os.path.join(output_folder, 'results.csv')
    savemodels=os.path.join(output_folder, 'models.pickle')
    
    # run
    _ = experiment_test_all_atn_predictors(dataset,
                                           target=target,
                                           covariates=covariates,
                                           stratify=stratify,
                                           repeats=repeats,
                                           splits=splits,
                                           seed=seed,
                                           savepath=savepath,
                                           savemodels=savemodels)
    
if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    main(**kwargs)