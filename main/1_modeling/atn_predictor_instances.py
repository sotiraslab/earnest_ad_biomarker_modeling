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
    BinaryManual('AmyloidComposite', 1.11, nickname='Amyloid Composite (SUVR>1.11)'),
    BinaryManual('Centiloid', 20, nickname='Centiloid (>20)'),
    BinaryGMM('AmyloidComposite', nickname='Amyloid Composite (GMM)'),
    BinaryGMM('Centiloid', nickname='Centiloid (GMM)'),
    BinaryZScore('AmyloidComposite', zcutoff=2.0, control_col='Control', nickname='Amyloid Composite (z>2.0)'),
    BinaryZScore('AmyloidComposite', zcutoff=2.5, control_col='Control', nickname='Amyloid Composite (z>2.5)'),
    BinaryZScore('Centiloid', zcutoff=2.0, control_col='Control', nickname='Centiloid (z>2.0)'),
    BinaryZScore('Centiloid', zcutoff=2.5, control_col='Control', nickname='Centiloid (z>2.5)')
    ]

AMYLOID_CATEGORICAL = [
    Quantiles('AmyloidComposite', nickname='Amyloid Composite (Quantiles)'),
    Quantiles('Centiloid', nickname='Centiloid (Quantiles)'),
    CategoricalStager(['MattssonEarlySUVR', 'MattssonIntermediateSUVR', 'MattssonLateSUVR'], nickname='Mattsson Staging'),
    CategoricalStager(collij_regions, groupings=collij_groupings, p=.5, nickname='Collij Staging')
    ]

AMYLOID_CONTINUOUS = [
    Continuous('AmyloidComposite', nickname='Amyloid Composite'),
    Continuous('Centiloid', nickname='Centiloid')
    ]

# ---- Tau (no-PVC) -----

TAU_BINARY = [
    BinaryGMM('META_TEMPORAL_TAU', nickname='MTT (GMM)'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.0, nickname='MTT (z>2.0)'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.5, nickname='MTT (z>2.5)'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.20, nickname='MTT (SUVR>1.20)'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.21, nickname='MTT (SUVR>1.21)'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.23, nickname='MTT (SUVR>1.23'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.33, nickname='MTT (SUVR>1.33)')
    ]

TAU_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAU', nickname='MTT (Quantiles)'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK34_TAU', 'BRAAK56_TAU'], nickname='Braak Staging (3)'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK3_TAU', 'BRAAK4_TAU', 'BRAAK5_TAU', 'BRAAK6_TAU'], nickname='Braak Staging (6)')
    ]

TAU_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAU', nickname='MTT'),
    Continuous('BRAAK1_TAU', nickname='Braak1'),
    Continuous('BRAAK34_TAU', nickname='Braak34'),
    Continuous('BRAAK56_TAU', nickname='Braak56')
    ]

# ---- GM Volume -----

GM_BINARY = [
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='Hippocampus (z<-2.0)'),
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='Hippocampus (z<-2.5)'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='MTV (z<-2.0)'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='MTV (z<-2.5)')
    ]

GM_CATEGORICAL = [
    Quantiles('HIPPOCAMPUS_VOL', nickname='Hippocampus (Quantiles)'),
    Quantiles('META_TEMPORAL_VOL', nickname='MTV (Quantiles)')
    ]

GM_CONTINUOUS = [
    Continuous('HIPPOCAMPUS_VOL', nickname='Hippocampus'),
    Continuous('META_TEMPORAL_VOL', nickname='MTV'),
    ]

# ---- Tau (PVC) ----

TAUPVC_BINARY = [
    BinaryGMM('META_TEMPORAL_TAUPVC', nickname='MTT (GMM) [PVC]'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.0, nickname='MTT (z>2.0) [PVC]'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.5, nickname='MTT (z>2.5) [PVC]'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.20, nickname='MTT (SUVR>1.20) [PVC]'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.21, nickname='MTT (SUVR>1.21) [PVC]'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.23, nickname='MTT (SUVR>1.23 [PVC]'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.33, nickname='MTT (SUVR>1.33) [PVC]')
    ]

TAUPVC_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAUPVC', nickname='MTT (Quantiles) [PVC]'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK34_TAUPVC', 'BRAAK56_TAUPVC'], nickname='Braak Staging (3) [PVC]'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK3_TAUPVC', 'BRAAK4_TAUPVC', 'BRAAK5_TAUPVC', 'BRAAK6_TAUPVC'], nickname='Braak Staging (6) [PVC]')
    ]

TAUPVC_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAUPVC', nickname='MTT [PVC]'),
    Continuous('BRAAK1_TAUPVC', nickname='Braak1 [PVC]'),
    Continuous('BRAAK34_TAUPVC', nickname='Braak34 [PVC]'),
    Continuous('BRAAK56_TAUPVC', nickname='Braak56 [PVC]')
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
    'gm': {
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