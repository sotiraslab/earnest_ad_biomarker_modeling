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

    features = [
        'DeltaPHC_MEM',
        'DeltaPHC_EXF',
        'DeltaPHC_VSP',
        'DeltaPHC_LAN',
        'DeltaPHC_PC1'
        ]

    for feature in features:
        # NOTE: Short has no effect
        # only included for consistency with other experiment scripts
        basename = os.path.basename(os.path.normpath(__file__))
        key = os.path.splitext(basename)[0].replace('NEUROPSYCH', feature)
        output_folder = setup_output(key)

        # EXPERIMENT PARAMETERS
        parameters = dict(
            dataset = pd.read_csv('datasets/maindata.csv'),
            target = feature,
            covariates=['Age', 'SexBinary', 'HasE4Binary'],
            stratify='CDRBinned',
            repeats=10,
            outer_splits=10,
            inner_splits=5,
            outer_seed=0,
            inner_seed=100,
            savepath=os.path.join(output_folder, 'results.csv'),
            savemodels=os.path.join(output_folder, 'models.pickle'),
            savelms=os.path.join(output_folder, 'lm.pickle')
            )

        _ = experiment_combo_atn_vs_baseline(**parameters)

if __name__ == '__main__':
    args = parse()
    kwargs = vars(args)
    main(**kwargs)
