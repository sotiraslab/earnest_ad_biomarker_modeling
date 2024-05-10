#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:08:53 2023

@author: earnestt1234
"""

from atn_modeling.atn_predictor_classes import (BinaryGMM,
                                                BinaryManual,
                                                BinaryZScore,
                                                CategoricalStager,
                                                GMMWithIndeterminateZone,
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
    BinaryManual('AmyloidComposite', 1.24, nickname='Su (SUVR>1.24)', atn='amyloid'),
    BinaryManual('AmyloidComposite', 1.42, nickname='Jack RW Amyloid', atn='amyloid'),
    BinaryManual('AmyloidComposite', 1.30, nickname='Jack Spec. Aymloid', atn='amyloid'),
    BinaryManual('Centiloid', 15, nickname='Centiloid (>15)', atn='amyloid'),
    BinaryManual('Centiloid', 20, nickname='Centiloid (>20)', atn='amyloid'),
    BinaryManual('Centiloid', 25, nickname='Centiloid (>25)', atn='amyloid'),
    BinaryManual('Centiloid', 30, nickname='Centiloid (>30)', atn='amyloid'),
    BinaryGMM('AmyloidComposite', nickname='Amyloid Composite (GMM)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.0, control_col='Control', nickname='Amyloid Composite (z>2.0)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.5, control_col='Control', nickname='Amyloid Composite (z>2.5)', atn='amyloid'),
    ]

AMYLOID_CATEGORICAL = [
    Quantiles('AmyloidComposite', nickname='Amyloid Composite (Quantiles)', atn='amyloid'),
    Quantiles('Centiloid', nickname='Centiloid (Quantiles)', atn='amyloid'),
    CategoricalStager(['MattssonEarlySUVR', 'MattssonIntermediateSUVR', 'MattssonLateSUVR'], nickname='Mattsson Staging', atn='amyloid'),
    CategoricalStager(collij_regions, groupings=collij_groupings, p=.5, nickname='Collij Staging', atn='amyloid'),
    GMMWithIndeterminateZone('AmyloidComposite', nickname='Amyloid Composite (UZ)', atn='amyloid'),
    GMMWithIndeterminateZone('Centiloid', nickname='Centiloid (UZ)', atn='amyloid'),
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
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.20, nickname='Jack Sens. MTT', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.21, nickname='Jack Spec. MTT', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.23, nickname='Jack Acc-Young MTT', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.33, nickname='Jack Acc-Matched MTT', atn='tau')
    ]

TAU_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAU', nickname='MTT (Quantiles)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK34_TAU', 'BRAAK56_TAU'], nickname='Braak Staging (3)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK3_TAU', 'BRAAK4_TAU', 'BRAAK5_TAU', 'BRAAK6_TAU'], nickname='Braak Staging (6)', atn='tau'),
    GMMWithIndeterminateZone('META_TEMPORAL_TAU', nickname='MTT (UZ)', atn='tau')
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
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.33, nickname='Jack Acc-Matched [PVC]', atn='taupvc')
    ]

TAUPVC_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAUPVC', nickname='MTT (Quantiles) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK34_TAUPVC', 'BRAAK56_TAUPVC'], nickname='Braak Staging (3) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK3_TAUPVC', 'BRAAK4_TAUPVC', 'BRAAK5_TAUPVC', 'BRAAK6_TAUPVC'], nickname='Braak Staging (6) [PVC]', atn='taupvc'),
    GMMWithIndeterminateZone('META_TEMPORAL_TAUPVC', nickname='MTT (UZ) [PVC]', atn='taupvc')
    ]

TAUPVC_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAUPVC', nickname='MTT [PVC]', atn='taupvc'),
    Continuous('BRAAK1_TAUPVC', nickname='Braak1 [PVC]', atn='taupvc'),
    Continuous('BRAAK34_TAUPVC', nickname='Braak34 [PVC]', atn='taupvc'),
    Continuous('BRAAK56_TAUPVC', nickname='Braak56 [PVC]', atn='taupvc')
    ]

# ---- CSF ----

CSF_AMYLOID_BINARY = [
    BinaryGMM('CSF_ABETA40', nickname='Aβ40 (GMM)', atn='csf_amyloid'),
    BinaryGMM('CSF_ABETA42', nickname='Aβ42 (GMM)', atn='csf_amyloid'),
    BinaryGMM('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (GMM)', atn='csf_amyloid'),
    ]

CSF_AMYLOID_CATEGORICAL = [
    GMMWithIndeterminateZone('CSF_ABETA40', nickname='Aβ40 (UZ)', atn='csf_amyloid'),
    GMMWithIndeterminateZone('CSF_ABETA42', nickname='Aβ42 (UZ)', atn='csf_amyloid'),
    GMMWithIndeterminateZone('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (UZ)', atn='csf_amyloid'),
    Quantiles('CSF_ABETA40', nickname='Aβ40 (Quantiles)', atn='csf_amyloid'),
    Quantiles('CSF_ABETA42', nickname='Aβ42 (Quantiles)', atn='csf_amyloid'),
    Quantiles('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (Quantiles)', atn='csf_amyloid'),
    ]

CSF_AMYLOID_CONTINUOUS = [
    Continuous('CSF_ABETA40', nickname='Aβ40', atn='csf_amyloid'),
    Continuous('CSF_ABETA42', nickname='Aβ42', atn='csf_amyloid'),
    Continuous('CSF_AB42OVER40', nickname='Aβ42/Aβ40', atn='csf_amyloid'),
    ]

CSF_TAU_BINARY = [
    BinaryGMM('CSF_TAU', nickname='tTau (GMM)', atn='csf_tau'),
    BinaryGMM('CSF_PTAU', nickname='pTau181 (GMM)', atn='csf_tau'),
    ]

CSF_TAU_CATEGORICAL = [
    GMMWithIndeterminateZone('CSF_TAU', nickname='tTau (UZ)', atn='csf_tau'),
    GMMWithIndeterminateZone('CSF_PTAU', nickname='pTau181 (UZ)', atn='csf_tau'),
    Quantiles('CSF_TAU', nickname='tTau (Quantiles', atn='csf_tau'),
    Quantiles('CSF_PTAU', nickname='pTau181 (Quantiles)', atn='csf_tau'),
    ]

CSF_TAU_CONTINUOUS = [
    Continuous('CSF_TAU', nickname='tTau', atn='csf_tau'),
    Continuous('CSF_PTAU', nickname='pTau181', atn='csf_tau'),
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

CSF_ATN_PREDICTORS = {
    'csf_amyloid': {
        'binary': CSF_AMYLOID_BINARY,
        'categorical': CSF_AMYLOID_CATEGORICAL,
        'continuous': CSF_AMYLOID_CONTINUOUS
        },
    'csf_tau': {
        'binary': CSF_TAU_BINARY,
        'categorical': CSF_TAU_CATEGORICAL,
        'continuous': CSF_TAU_CONTINUOUS
        },
    }

ATN_PREDICTORS_PLUS_CSF = ATN_PREDICTORS.copy()
ATN_PREDICTORS_PLUS_CSF.update(CSF_ATN_PREDICTORS)

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

CSF_ATN_PREDICTORS_FLAT = [
    *CSF_AMYLOID_BINARY,
    *CSF_AMYLOID_CATEGORICAL,
    *CSF_AMYLOID_CONTINUOUS,
    *CSF_TAU_BINARY,
    *CSF_TAU_CATEGORICAL,
    *CSF_TAU_CONTINUOUS,]

# Accessor functions -------

def get_models_by_nickname(nicknames):
    if isinstance(nicknames, str):
        nicknames = [nicknames]
    return [m for m in ATN_PREDICTORS_FLAT if m.nickname in nicknames]

# Check for duplicate names -------

nicknames = [m.nickname for m in ATN_PREDICTORS_FLAT] + [m.nickname for m in CSF_ATN_PREDICTORS_FLAT]
nicknames_set = set(nicknames)
assert len(nicknames) == len(nicknames_set), 'Resolve duplicated nicknames for ATN predictors!'
