#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:12:35 2024

@author: tom.earnest
"""
import os

import pandas as pd

from atn_modeling.experiments import experiment_combo_atn_vs_baseline_PLASMA
from common import parse, setup_output

def main(short=False):

    # NOTE: Short has no effect
    # only included for consistency with other experiment scripts

    output_folder = setup_output(__file__)

    # parameters
    dataset=pd.read_csv('datasets/maindata_plasma.csv')
    target='DeltaMMSE'
    covariates=['Age', 'SexBinary', 'HasE4Binary']
    stratify='CDRBinned'
    repeats=10
    outer_splits=10
    inner_splits=5
    outer_seed=0
    inner_seed=100
    savepath=os.path.join(output_folder, 'results.csv')
    savemodels=os.path.join(output_folder, 'models.pickle')
    savelms=os.path.join(output_folder, 'lm.pickle')
    testing_fiilter=lambda x: x.copy().loc[x['CDRBinned'].eq('0.0')]

    # run
    _ = experiment_combo_atn_vs_baseline_PLASMA(dataset,
                                                target=target,
                                                covariates=covariates,
                                                stratify=stratify,
                                                repeats=repeats,
                                                outer_splits=outer_splits,
                                                inner_splits=inner_splits,
                                                outer_seed=outer_seed,
                                                inner_seed=inner_seed,
                                                savepath=savepath,
                                                savemodels=savemodels,
                                                savelms=savelms,
                                                testing_filter=testing_fiilter)

if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    main(**kwargs)
