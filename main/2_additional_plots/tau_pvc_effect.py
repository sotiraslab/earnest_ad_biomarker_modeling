#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 09:51:00 2023

@author: tom.earnest
"""

# ---- imports
import os
import pickle

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

import sys
sys.path.append('../1_modeling')
# from helpers impo

# ---- load saved models from binary combination experiment

import sys
sys.path.append('../1_modeling')

separate_models_path = '../../outputs/exp0_individual_atn_models_global_cognition/models.pickle'
combo_models_path = '../../outputs/exp2_combo_atn_models_global_cognition/models.pickle'
with open(separate_models_path, 'rb') as f:
    sep_models = pickle.load(f)
    
with open(combo_models_path, 'rb') as f:
    com_models = pickle.load(f)
    
# ---- select all tau models

pvc_model_names = []
nopvc_model_names = []

for k, v in sep_models.items():
    m = v[0]
    if m.atn == 'tau':
        nopvc_model_names.append(k)
    elif m.atn == 'taupvc':
        pvc_model_names.append(k)
    else:
        pass

#%%
