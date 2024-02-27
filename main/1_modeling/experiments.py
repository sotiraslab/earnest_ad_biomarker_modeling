#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 15:16:00 2024

@author: tom.earnest
"""

from atn_predictor_instances import (AMYLOID_BINARY,
                                     AMYLOID_CATEGORICAL,
                                     AMYLOID_CONTINUOUS,
                                     TAU_BINARY,
                                     TAU_CATEGORICAL,
                                     TAU_CONTINUOUS,
                                     GM_BINARY,
                                     GM_CATEGORICAL,
                                     GM_CONTINUOUS)

def experiment_test_all_predictors_LM(dataset, target,
                                      covariates=['Age', 'SexBinary', 'HasE4Binary'],
                                      statify='CDRBinned',
                                      splits=5,
                                      repeats=20,
                                      seed=0):
    pass