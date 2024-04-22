#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import os

import numpy as np
import pandas as pd

from atn_modeling.experiments import experiment_svm
from common import parse, setup_output

def main(short=False):
    
    # setup output
    output_folder = setup_output(__file__, short=short)

    # EXPERIMENT PARAMETERS
    params = dict(
        dataset=pd.read_csv('datasets/maindata_long.csv'),
        target='DeltaADSP',
        covariates=['Age', 'SexBinary', 'HasE4Binary'],
        stratify='CDRBinned',
        search_C_linear = [2 ** -7] if short else list(2. ** np.arange(-10, 1, 1)),
        search_C_rbf = [2] if short else list(2. ** np.arange(-5, 17, 2)),
        search_gamma = [2 ** -11] if short else list(2. ** np.arange(-15, 5, 2)),
        search_kernel = ['linear', 'rbf'],
        repeats=10,
        outer_splits=10,
        inner_splits=5,
        outer_seed=0,
        inner_seed=100,
        savepath=os.path.join(output_folder, 'results.csv'),
        savemodels=os.path.join(output_folder, 'models.pickle')
        )
    
    # run
    _ = experiment_svm(**params)
    
if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    main(**kwargs)