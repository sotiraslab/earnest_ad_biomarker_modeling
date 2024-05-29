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
    BinaryManual('AmyloidComposite', 1.11, nickname='Aβ SUVR>1.11', atn='amyloid'),
    BinaryManual('AmyloidComposite', 1.24, nickname='Aβ SUVR>1.24', atn='amyloid'),
    BinaryManual('AmyloidComposite', 1.42, nickname='Aβ SUVR>1.42', atn='amyloid'),
    BinaryManual('AmyloidComposite', 1.30, nickname='Aβ SUVR>1.30', atn='amyloid'),
    BinaryManual('Centiloid', 15, nickname='Centiloid>15', atn='amyloid'),
    BinaryManual('Centiloid', 20, nickname='Centiloid>20', atn='amyloid'),
    BinaryManual('Centiloid', 25, nickname='Centiloid>25', atn='amyloid'),
    BinaryManual('Centiloid', 30, nickname='Centiloid>30', atn='amyloid'),
    BinaryGMM('AmyloidComposite', nickname='Aβ Composite (GMM)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.0, control_col='Control', nickname='Aβ Composite (z>2.0)', atn='amyloid'),
    BinaryZScore('AmyloidComposite', zcutoff=2.5, control_col='Control', nickname='Aβ Composite (z>2.5)', atn='amyloid'),
    ]

AMYLOID_CATEGORICAL = [
    Quantiles('AmyloidComposite', nickname='Aβ Composite (Quartiles)', atn='amyloid'),
    Quantiles('Centiloid', nickname='Centiloid (Quartiles)', atn='amyloid'),
    CategoricalStager(['MattssonEarlySUVR', 'MattssonIntermediateSUVR', 'MattssonLateSUVR'], nickname='Mattsson Staging', atn='amyloid'),
    CategoricalStager(collij_regions, groupings=collij_groupings, p=.5, nickname='Collij Staging', atn='amyloid'),
    GMMWithIndeterminateZone('AmyloidComposite', nickname='Aβ Composite (BIZ)', atn='amyloid'),
    GMMWithIndeterminateZone('Centiloid', nickname='Centiloid (BIZ)', atn='amyloid'),
    ]

AMYLOID_CONTINUOUS = [
    Continuous('AmyloidComposite', nickname='Aβ Composite', atn='amyloid'),
    Continuous('Centiloid', nickname='Centiloid', atn='amyloid')
    ]

# ---- Tau (no-PVC) -----

TAU_BINARY = [
    BinaryGMM('META_TEMPORAL_TAU', nickname='MT tau (GMM)', atn='tau'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.0, nickname='MT tau (z>2.0)', atn='tau'),
    BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.5, nickname='MT tau (z>2.5)', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.20, nickname='Tau SUVR>1.20', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.21, nickname='Tau SUVR>1.21', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.23, nickname='Tau SUVR>1.23', atn='tau'),
    BinaryManual('META_TEMPORAL_TAU', cutoff=1.33, nickname='Tau SUVR>1.33', atn='tau')
    ]

TAU_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAU', nickname='MT tau (Quartiles)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK34_TAU', 'BRAAK56_TAU'], nickname='Braak staging (3)', atn='tau'),
    CategoricalStager(['BRAAK1_TAU', 'BRAAK3_TAU', 'BRAAK4_TAU', 'BRAAK5_TAU', 'BRAAK6_TAU'], nickname='Braak staging (6)', atn='tau'),
    GMMWithIndeterminateZone('META_TEMPORAL_TAU', nickname='MT tau (BIZ)', atn='tau')
    ]

TAU_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAU', nickname='MT tau SUVR', atn='tau'),
    Continuous('BRAAK1_TAU', nickname='Braak I SUVR', atn='tau'),
    Continuous('BRAAK34_TAU', nickname='Braak III/IV SUVR', atn='tau'),
    Continuous('BRAAK56_TAU', nickname='Braak V/VI SUVR', atn='tau')
    ]

# ---- GM Volume -----

GM_BINARY = [
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='Hippocampus (z<-2.0)', atn='neurodegeneration'),
    BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='Hippocampus (z<-2.5)', atn='neurodegeneration'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='MT volume (z<-2.0)', atn='neurodegeneration'),
    BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='MT volume (z<-2.5)', atn='neurodegeneration')
    ]

GM_CATEGORICAL = [
    Quantiles('HIPPOCAMPUS_VOL', nickname='Hippocampus (Quartiles)', atn='neurodegeneration'),
    Quantiles('META_TEMPORAL_VOL', nickname='MTV (Quartiles)', atn='neurodegeneration')
    ]

GM_CONTINUOUS = [
    Continuous('HIPPOCAMPUS_VOL', nickname='Hippocampus', atn='neurodegeneration'),
    Continuous('META_TEMPORAL_VOL', nickname='MT volume', atn='neurodegeneration'),
    ]

# ---- Tau (PVC) ----

TAUPVC_BINARY = [
    BinaryGMM('META_TEMPORAL_TAUPVC', nickname='MT tau (GMM) [PVC]', atn='taupvc'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.0, nickname='MT tau (z>2.0) [PVC]', atn='taupvc'),
    BinaryZScore('META_TEMPORAL_TAUPVC', 'Control', zcutoff=2.5, nickname='MT tau (z>2.5) [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.20, nickname='Tau SUVR>1.20 [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.21, nickname='Tau SUVR>1.21 [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.23, nickname='Tau SUVR>1.23 [PVC]', atn='taupvc'),
    BinaryManual('META_TEMPORAL_TAUPVC', cutoff=1.33, nickname='Tau SUVR>1.33 [PVC]', atn='taupvc')
    ]

TAUPVC_CATEGORICAL = [
    Quantiles('META_TEMPORAL_TAUPVC', nickname='MT tau (Quartiles) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK34_TAUPVC', 'BRAAK56_TAUPVC'], nickname='Braak staging (3) [PVC]', atn='taupvc'),
    CategoricalStager(['BRAAK1_TAUPVC', 'BRAAK3_TAUPVC', 'BRAAK4_TAUPVC', 'BRAAK5_TAUPVC', 'BRAAK6_TAUPVC'], nickname='Braak staging (6) [PVC]', atn='taupvc'),
    GMMWithIndeterminateZone('META_TEMPORAL_TAUPVC', nickname='MT tau (BIZ) [PVC]', atn='taupvc')
    ]

TAUPVC_CONTINUOUS = [
    Continuous('META_TEMPORAL_TAUPVC', nickname='MT tau [PVC]', atn='taupvc'),
    Continuous('BRAAK1_TAUPVC', nickname='Braak I SUVR [PVC]', atn='taupvc'),
    Continuous('BRAAK34_TAUPVC', nickname='Braak III/IV SUVR [PVC]', atn='taupvc'),
    Continuous('BRAAK56_TAUPVC', nickname='Braak V/VI [PVC]', atn='taupvc')
    ]

# ---- CSF ----

CSF_AMYLOID_BINARY = [
    BinaryGMM('CSF_ABETA40', nickname='Aβ40 (GMM)', atn='csf_amyloid'),
    BinaryGMM('CSF_ABETA42', nickname='Aβ42 (GMM)', atn='csf_amyloid'),
    BinaryGMM('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (GMM)', atn='csf_amyloid'),
    ]

CSF_AMYLOID_CATEGORICAL = [
    GMMWithIndeterminateZone('CSF_ABETA40', nickname='Aβ40 (BIZ)', atn='csf_amyloid'),
    GMMWithIndeterminateZone('CSF_ABETA42', nickname='Aβ42 (BIZ)', atn='csf_amyloid'),
    GMMWithIndeterminateZone('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (BIZ)', atn='csf_amyloid'),
    Quantiles('CSF_ABETA40', nickname='Aβ40 (Quartiles)', atn='csf_amyloid'),
    Quantiles('CSF_ABETA42', nickname='Aβ42 (Quartiles)', atn='csf_amyloid'),
    Quantiles('CSF_AB42OVER40', nickname='Aβ42/Aβ40 (Quartiles)', atn='csf_amyloid'),
    ]

CSF_AMYLOID_CONTINUOUS = [
    Continuous('CSF_ABETA40', nickname='Aβ40', atn='csf_amyloid'),
    Continuous('CSF_ABETA42', nickname='Aβ42', atn='csf_amyloid'),
    Continuous('CSF_AB42OVER40', nickname='Aβ42/Aβ40', atn='csf_amyloid'),
    ]

CSF_TAU_BINARY = [
    BinaryGMM('CSF_PTAU', nickname='pTau181 (GMM)', atn='csf_tau'),
    ]

CSF_TAU_CATEGORICAL = [
    GMMWithIndeterminateZone('CSF_PTAU', nickname='pTau181 (BIZ)', atn='csf_tau'),
    Quantiles('CSF_PTAU', nickname='pTau181 (Quartiles)', atn='csf_tau'),
    ]

CSF_TAU_CONTINUOUS = [
    Continuous('CSF_PTAU', nickname='pTau181', atn='csf_tau'),
    ]

CSF_NEURODEGENERATION_BINARY = [
    BinaryGMM('CSF_TAU', nickname='tTau (GMM)', atn='csf_neurodegeneration'),
    ]

CSF_NEURODEGENERATION_CATEGORICAL = [
    GMMWithIndeterminateZone('CSF_TAU', nickname='tTau (BIZ)', atn='csf_neurodegeneration'),
    Quantiles('CSF_TAU', nickname='tTau (Quartiles)', atn='csf_neurodegeneration'),
    ]

CSF_NEURODEGENERATION_CONTINUOUS = [
    Continuous('CSF_TAU', nickname='tTau', atn='csf_neurodegeneration')
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
    'csf_neurodegeneration': {
        'binary': CSF_NEURODEGENERATION_BINARY,
        'categorical': CSF_NEURODEGENERATION_CATEGORICAL,
        'continuous': CSF_NEURODEGENERATION_CONTINUOUS
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
    *CSF_TAU_CONTINUOUS,
    *CSF_NEURODEGENERATION_BINARY,
    *CSF_NEURODEGENERATION_CATEGORICAL,
    *CSF_NEURODEGENERATION_CONTINUOUS,]

ALL_PREDICTORS_FLAT = ATN_PREDICTORS_FLAT + CSF_ATN_PREDICTORS_FLAT

# Accessor functions -------

def get_models_by_nickname(nicknames):
    if isinstance(nicknames, str):
        nicknames = [nicknames]
    return [m for m in ALL_PREDICTORS_FLAT if m.nickname in nicknames]

# Check for duplicate names -------

nicknames = [m.nickname for m in ALL_PREDICTORS_FLAT]
nicknames_set = set(nicknames)
assert len(nicknames) == len(nicknames_set), 'Resolve duplicated nicknames for ATN predictors!'
