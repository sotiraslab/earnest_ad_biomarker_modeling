#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""

import os

import pandas as pd

from atn_modeling.experiments import experiment_combo_atn_vs_baseline
from common import parse, setup_output

def main(short=False):
    
    # NOTE: Short has no effect
    # only included for consistency with other experiment scripts

    # setup output
    output_folder = setup_output(__file__)

    # parameters
    params = dict(
        dataset=pd.read_csv('datasets/maindata.csv'),
        target='PHC_GLOBAL',
        covariates=['Age', 'SexBinary', 'HasE4Binary'],
        stratify='CDRBinned',
        repeats=10,
        outer_splits=10,
        inner_splits=5,
        outer_seed=0,
        inner_seed=100,
        savepath=os.path.join(output_folder, 'results.csv'),
        savemodels=os.path.join(output_folder, 'models.pickle'),
        savelms=os.path.join(output_folder, 'lm.pickle'),
        coefficient_analysis=True
        )

    # run
    _ = experiment_combo_atn_vs_baseline(**params)

if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    main(**kwargs)
