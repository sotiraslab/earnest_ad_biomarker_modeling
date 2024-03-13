#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:08:53 2023

@author: earnestt1234
"""

from atn_predictor_classes import (BinaryGMM,
                                   BinaryManual,
                                   BinaryZScore,
                                   CategoricalStager,
                                   Continuous,
                                   Quantiles)

collij_regions = [
    "AV45_CTX_TOT_POSTERIORCINGULATE_SUVR",
    "AV45_CTX_TOT_ISTHMUSCINGULATE_SUVR",
    "CollijAnteriorCingulate",
    
    "AV45_CTX_TOT_LATERALORBITOFRONTAL_SUVR",
    "AV45_CTX_TOT_PARACENTRAL_SUVR",
    "AV45_CTX_TOT_PRECUNEUS_SUVR",
    "AV45_CTX_TOT_MEDIALORBITOFRONTAL_SUVR",
    
    "AV45_CTX_TOT_INSULA_SUVR",
    "AV45_CTX_TOT_FUSIFORM_SUVR",
    "AV45_CTX_TOT_PRECENTRAL_SUVR",
    "AV45_CTX_TOT_INFERIORTEMPORAL_SUVR",
    "AV45_CTX_TOT_PARAHIPPOCAMPAL_SUVR",
    "CollijInferiorFrontal",
    "AV45_CTX_TOT_SUPERIORFRONTAL_SUVR",
    "AV45_CTX_TOT_LINGUAL_SUVR",
    "AV45_CTX_TOT_SUPRAMARGINAL_SUVR",
    "AV45_CTX_TOT_INFERIORPARIETAL_SUVR",
    "AV45_CTX_TOT_CUNEUS_SUVR",
    "CollijMiddleFrontal",
    
    "AV45_CTX_TOT_LATERALOCCIPITAL_SUVR",
    "AV45_CTX_TOT_SUPERIORPARIETAL_SUVR",
    "AV45_CTX_TOT_MIDDLETEMPORAL_SUVR",
    "AV45_CTX_TOT_SUPERIORTEMPORAL_SUVR",
    "AV45_CTX_TOT_POSTCENTRAL_SUVR",
    "AV45_CTX_TOT_ENTORHINAL_SUVR",
    "AV45_CTX_TOT_FRONTALPOLE_SUVR",
    "AV45_CTX_TOT_TEMPORALPOLE_SUVR"
    ]

collij_groupings = [0] * 3 + [1] * 4 + [2] * 12 + [3] * 8


# ---- Amyloid -----
AMYLOID_BINARY = [
    BinaryManual('AmyloidComposite', 1.11, nickname='Landau (SUVR>1.11)', atn='amyloid'),
    BinaryManual('Centiloid', 20, nickname='Centiloid (>20)', atn='amyloid'),
    BinaryGMM('AmyloidComposite', nickname='Amyloid Composite (GMM)', atn='amyloid'),
    BinaryGMM('Centiloid', nickname='Centiloid (GMM)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.0, control_col='Control', nickname='Amyloid Composite (z>2.0)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.5, control_col='Control', nickname='Amyloid Composite (z>2.5)', atn='amyloid'),
    BinaryZScore('Centiloid', zcutoff=2.0, control_col='Control', nickname='Centiloid (z>2.0)', atn='amyloid'),
    BinaryZScore('Centiloid', zcutoff=2.5, control_col='Control', nickname='Centiloid (z>2.5)', atn='amyloid')
    ]

AMYLOID_CATEGORICAL = [
    Quantiles('AmyloidComposite', nickname='Amyloid Composite (Quantiles)', atn='amyloid'),
    Quantiles('Centiloid', nickname='Centiloid (Quantiles)', atn='amyloid'),
    CategoricalStager(['MattssonEarlySUVR', 'MattssonIntermediateSUVR', 'MattssonLateSUVR'], nickname='Mattsson Staging', atn='amyloid'),
    CategoricalStager(collij_regions, groupings=collij_groupings, p=.5, nickname='Collij Staging', atn='amyloid')
    ]

AMYLOID_CONTINUOUS = [
    Continuous('AmyloidComposite', nickname='Amyloid Composite', atn='amyloid'),
    Continuous('Centiloid', nickname='Centiloid', atn='amyloid')
    ]

# ---- Tau (no-PVC) -----

TAU_BINARY = [
    BinaryGMM('META_TEMPORAL_TAU', nickname='MTT (GMM)', atn='tau'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.0, nickname='MTT (z>2.0)', atn='tau'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.5, nickname='MTT (z>2.5)', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.20, nickname='Jack Sens.', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.21, nickname='Jack Spec.', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.23, nickname='Jack Acc-Young', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.33, nickname='Jack Acc-Matched', atn='tau')
    ]

TAU_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAU', nickname='MTT (Quantiles)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK34_TAU', 'BRAAK56_TAU'], nickname='Braak Staging (3)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK3_TAU', 'BRAAK4_TAU', 'BRAAK5_TAU', 'BRAAK6_TAU'], nickname='Braak Staging (6)', atn='tau')
    ]

TAU_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAU', nickname='MTT', atn='tau'),
    Continuous('BRAAK1_TAU', nickname='Braak1', atn='tau'),
    Continuous('BRAAK34_TAU', nickname='Braak34', atn='tau'),
    Continuous('BRAAK56_TAU', nickname='Braak56', atn='tau')
    ]

# ---- GM Volume -----

GM_BINARY = [
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='Hippocampus (z<-2.0)', atn='neurodegeneration'),
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='Hippocampus (z<-2.5)', atn='neurodegeneration'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='MTV (z<-2.0)', atn='neurodegeneration'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='MTV (z<-2.5)', atn='neurodegeneration')
    ]

GM_CATEGORICAL = [
    Quantiles('HIPPOCAMPUS_VOL', nickname='Hippocampus (Quantiles)', atn='neurodegeneration'),
    Quantiles('META_TEMPORAL_VOL', nickname='MTV (Quantiles)', atn='neurodegeneration')
    ]

GM_CONTINUOUS = [
    Continuous('HIPPOCAMPUS_VOL', nickname='Hippocampus', atn='neurodegeneration'),
    Continuous('META_TEMPORAL_VOL', nickname='MTV', atn='neurodegeneration'),
    ]

# ---- Tau (PVC) ----

TAUPVC_BINARY = [
    BinaryGMM('META_TEMPORAL_TAUPVC', nickname='MTT (GMM) [PVC]', atn='taupvc'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.0, nickname='MTT (z>2.0) [PVC]', atn='taupvc'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.5, nickname='MTT (z>2.5) [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.20, nickname='Jack Sens. [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.21, nickname='Jack Spec. [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.23, nickname='Jack Acc-Young [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.33, nickname='MJack Acc-Matched [PVC]', atn='taupvc')
    ]

TAUPVC_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAUPVC', nickname='MTT (Quantiles) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK34_TAUPVC', 'BRAAK56_TAUPVC'], nickname='Braak Staging (3) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK3_TAUPVC', 'BRAAK4_TAUPVC', 'BRAAK5_TAUPVC', 'BRAAK6_TAUPVC'], nickname='Braak Staging (6) [PVC]', atn='taupvc')
    ]

TAUPVC_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAUPVC', nickname='MTT [PVC]', atn='taupvc'),
    Continuous('BRAAK1_TAUPVC', nickname='Braak1 [PVC]', atn='taupvc'),
    Continuous('BRAAK34_TAUPVC', nickname='Braak34 [PVC]', atn='taupvc'),
    Continuous('BRAAK56_TAUPVC', nickname='Braak56 [PVC]', atn='taupvc')
    ]

# ---- Big dictionary -----

ATN_PREDICTORS = {
    'amyloid': {
        'binary': AMYLOID_BINARY,
        'categorical': AMYLOID_CATEGORICAL,
        'continuous': AMYLOID_CONTINUOUS
        },
    'tau': {
        'binary': TAU_BINARY,
        'categorical': TAU_CATEGORICAL,
        'continuous': TAU_CONTINUOUS
        }, 
    'neurodegeneration': {
        'binary': GM_BINARY,
        'categorical': GM_CATEGORICAL,
        'continuous': GM_CONTINUOUS
        },
    'taupvc': {
        'binary': TAUPVC_BINARY,
        'categorical': TAUPVC_CATEGORICAL,
        'continuous': TAUPVC_CONTINUOUS
        }, 
    }

ATN_PREDICTORS_FLAT = [*AMYLOID_BINARY,
                       *AMYLOID_CATEGORICAL,
                       *AMYLOID_CONTINUOUS,
                       *TAU_BINARY,
                       *TAU_CATEGORICAL,
                       *TAU_CONTINUOUS,
                       *GM_BINARY,
                       *GM_CATEGORICAL,
                       *GM_CONTINUOUS,
                       *TAUPVC_BINARY,
                       *TAUPVC_CATEGORICAL,
                       *TAUPVC_CONTINUOUS]

# Accessor functions -------

def get_models_by_nickname(nicknames):
    if isinstance(nicknames, str):
        nicknames = [nicknames]
    return [m for m in ATN_PREDICTORS_FLAT if m.nickname in nicknames]