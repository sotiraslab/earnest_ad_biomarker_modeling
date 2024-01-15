#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:08:53 2023

@author: earnestt1234
"""

import inspect
import operator

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import root_scalar
from sklearn.mixture import GaussianMixture
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR, LinearSVR

# Base Class -----

class ATNPredictor:
    
    def __init__(self, nickname=None):
        self.nickname = nickname

    def __repr__(self):
        names = inspect.getfullargspec(self.__init__).args
        d = {}
        for name in names:
            try:
                value = getattr(self, name)
                d[name] = value
            except AttributeError:
                continue
        kv = ', '.join([f'{k}="{v}"' if isinstance(v, str) else f'{k}={v}'
                        for k, v in d.items()])
        classname = self.__class__.__name__
        s = f"{classname}({kv})"

        return s
    
    @property
    def name(self):
        return self.nickname if self.nickname else repr(self)

# Binary predictors ----

class BinaryZScore(ATNPredictor):

    def __init__(self, y_col, control_col, zcutoff=2.0, greater=True, nickname=None):
        super().__init__(nickname)
        self.y_col = y_col
        self.control_col = control_col
        self.zcutoff = zcutoff
        self.operator = operator.ge if greater else operator.le

        self.cutoff = None
        self.mean = None
        self.std = None

    def fit(self, data):
        controls = data.loc[data[self.control_col].astype(bool),  self.y_col]
        self.mean = controls.mean()
        self.std = controls.std()
        self.cutoff = (self.zcutoff * self.std) + self.mean

    def covariates(self, data):
        y = data[self.y_col]
        z = (y - self.mean) / self.std
        return np.where(self.operator(z, self.zcutoff), 1., 0.)[:, np.newaxis]

class BinaryGMM(ATNPredictor):

    def __init__(self, y_col, greater=True, nickname=None):
        super().__init__(nickname)
        self.y_col = y_col
        self.operator = operator.ge if greater else operator.le

        self.gmm = None
        self.cutoff = None

    def fit(self, data):
        self.gmm = GaussianMixture(n_components=2)
        self.gmm.fit(data[self.y_col].to_numpy()[:, np.newaxis])

        means = self.gmm.means_
        stds = np.sqrt(self.gmm.covariances_.flatten())
        weights = self.gmm.weights_

        g1 = lambda x: norm.pdf(x, loc = means[0], scale = stds[0]) * weights[0]
        g2 = lambda x: norm.pdf(x, loc = means[1], scale = stds[1]) * weights[1]
        diff = lambda x: g2(x) - g1(x)
        try:
            self.cutoff = root_scalar(diff, bracket=list(means)).root
        except ValueError:
            self.cutoff = (means[0] + 2*stds[0])[0]

    def covariates(self, data):
        y = data[self.y_col]
        return np.where(self.operator(y, self.cutoff), 1., 0.)[:, np.newaxis]

class BinaryManual(ATNPredictor):

    def __init__(self, y_col, cutoff, greater=True, nickname=None):
        super().__init__(nickname)
        self.y_col = y_col
        self.cutoff = cutoff
        self.operator = operator.ge if greater else operator.le

    def fit(self, data):
        pass

    def covariates(self, data):
        y = data[self.y_col]
        return np.where(self.operator(y, self.cutoff), 1., 0.)[:, np.newaxis]

# Categorical predictors ----

def assign_frequency_stage(data, groupings=None, p='any', atypical='NS'):

    if groupings == None:
        groupings = list(range(data.shape[1]))

    unique_stages = sorted(list(set(groupings)))
    n = len(unique_stages)
    stage_mat = np.zeros((len(data), n))

    for i in unique_stages:
        sub = data[:, np.array(groupings) == i]
        freqs = sub.sum(axis=1) / data.shape[1]

        if p == 'any':
            positive = freqs > 0
        elif p == 'all':
            positive = freqs == 1
        else:
            positive = freqs >= p

        stage_mat[:, i] = positive

    diffs = np.diff(stage_mat, axis=1)
    if n == 2:
        increasing = diffs <= 0
    else:
        increasing = np.all(diffs <= 0, axis=1)
    stage = np.where(increasing, stage_mat.sum(axis=1).astype(int), atypical)

    cats = [str(i) for i in range(0, len(unique_stages) + 1)] + [str(atypical)]
    return pd.Categorical(stage, categories=cats)

class Quantiles(ATNPredictor):

    def __init__(self, y_col, nickname=None):
        super().__init__(nickname)
        self.y_col = y_col

        self.quantiles = None

    def fit(self, data):
        self.quantiles = np.quantile(data[self.y_col], [.25, .5, .75, 1])

    def covariates(self, data):
        return np.digitize(data[self.y_col], self.quantiles)[:, np.newaxis]

class CategoricalStager(ATNPredictor):

    def __init__(self, columns, groupings=None, method='gmm', non_stageable='NS',
                 p='any', nickname=None, **kwargs):
        super().__init__(nickname)
        self.columns = columns
        self.groupings = groupings
        self.method = method
        self.non_stageable = non_stageable
        self.p = p
        self.kwargs = kwargs

        self.binary_models = None

    def fit(self, data):

        if self.method == 'gmm':
            model_class = BinaryGMM
        else:
            raise NotImplementedError(f'method "{self.method}" not implemented')

        self.binary_models = []

        for col in self.columns:
            model = model_class(y_col=col, **self.kwargs)
            model.fit(data)
            self.binary_models.append(model)

    def covariates(self, data):

        binary = np.zeros((len(data), len(self.columns)))
        for i, col in enumerate(self.binary_models):
            m = self.binary_models[i]
            binary[:, i] = m.covariates(data)[:, 0]

        stages = assign_frequency_stage(binary, groupings=self.groupings,
                                        p=self.p, atypical=self.non_stageable)
        return pd.get_dummies(stages).to_numpy()

# Continuous predictors ----

class Continuous(ATNPredictor):

    def __init__(self, y_col, nickname=None):
        super().__init__(nickname)
        self.y_col = y_col

    def fit(self, data):
        pass

    def covariates(self, data):
        return data[self.y_col].to_numpy()[:, np.newaxis]

# Multivariate -----

class MultivariateSVR:

    def __init__(self, predictors, target, **kwargs):
        self.predictors = predictors
        self.target = target

        kernel = kwargs.get('kernel', None)
        if kernel and kernel == 'linear':
            del kwargs['kernel']
            SVM = LinearSVR(max_iter=int(1e5), **kwargs)
        else:
            SVM = SVR(**kwargs)
        self.kwargs = kwargs

        self.pipeline = Pipeline([('scaler', StandardScaler()),
                                  ('svm', SVM)])

    def fit(self, data):
        X = data[self.predictors].to_numpy()
        y = data[self.target].to_numpy()

        self.pipeline.fit(X, y)

    def predict(self, data):
        X = data[self.predictors].to_numpy()
        return self.pipeline.predict(X)

# ATN predictory dictionary -----

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

ATN_PREDICTORS_DICT = {
    'amyloid': {
        'binary': {
            'composite_1.11': BinaryManual('AMYLOID_COMPOSITE', 1.11, nickname='Amyloid Composite (SUVR>1.11)'),
            'centiloid_20': BinaryManual('Centiloid', 20, nickname='Centiloid (>20)'),
            'composite_gmm': BinaryGMM('AMYLOID_COMPOSITE', nickname='Amyloid Composite (GMM)'),
            'centiloid_gmm': BinaryGMM('Centiloid', nickname='Centiloid (GMM)'),
            'composite_z2.0': BinaryZScore('AMYLOID_COMPOSITE', zcutoff=2.0, control_col='Control', nickname='Amyloid Composite (z>2.0)'),
            'composite_z2.5': BinaryZScore('AMYLOID_COMPOSITE', zcutoff=2.5, control_col='Control', nickname='Amyloid Composite (z>2.5)'),
            'centiloid_z2.0': BinaryZScore('Centiloid', zcutoff=2.0, control_col='Control', nickname='Centiloid (z>2.0)'),
            'centiloid_z2.5': BinaryZScore('Centiloid', zcutoff=2.5, control_col='Control', nickname='Centiloid (z>2.5)')},
        'categorical': {
            'composite_quantiles': Quantiles('AMYLOID_COMPOSITE', nickname='Amyloid Composite (Quantiles)'),
            'centiloid_quantiles': Quantiles('Centiloid', nickname='Centiloid (Quantiles)',),
            'mattsson_staging_gmm': CategoricalStager(['MattssonEarlySUVR', 'MattssonIntermediateSUVR', 'MattssonLateSUVR'], nickname='Mattsson Staging'),
            'collij_staging_gmm': CategoricalStager(collij_regions, groupings=collij_groupings, p=.5, nickname='Collij Staging')},
        'continuous': {
            'composite': Continuous('AMYLOID_COMPOSITE', nickname='Amyloid Composite (Quantiles)'),
            'centiloid': Continuous('Centiloid', nickname='Centiloid (Quantiles)')}},
    'tau': {
        'binary': {
            'mtt_gmm': BinaryGMM('META_TEMPORAL_TAU', nickname='MTT (GMM)'),
            'mtt_z2.0': BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.0, nickname='MTT (z>2.0)'),
            'mtt_z2.5': BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.5, nickname='MTT (z>2.5)'),
            'mtt_1.20': BinaryManual('META_TEMPORAL_TAU', cutoff=1.20, nickname='MTT (SUVR>1.20)'),
            'mtt_1.21': BinaryManual('META_TEMPORAL_TAU', cutoff=1.21, nickname='MTT (SUVR>1.21)'),
            'mtt_1.23': BinaryManual('META_TEMPORAL_TAU', cutoff=1.23, nickname='MTT (SUVR>1.23'),
            'mtt_1.33':  BinaryManual('META_TEMPORAL_TAU', cutoff=1.33, nickname='MTT (SUVR>1.33)')},
        'categorical': {
            'mtt_quantiles': Quantiles('META_TEMPORAL_TAU', nickname='MTT (Quantiles)'),
            'braak_stage_gmm': CategoricalStager(['BRAAK1_TAU', 'BRAAK34_TAU', 'BRAAK56_TAU'], nickname='Braak Staging')},
        'continuous': {
            'mtt': Continuous('META_TEMPORAL_TAU', nickname='MTT'),
            'braak1': Continuous('BRAAK1_TAU', nickname='Braak1'),
            'braak34': Continuous('BRAAK34_TAU', nickname='Braak34'),
            'braak56': Continuous('BRAAK56_TAU', nickname='Braak56')}
        },
    'neurodegeneration': {
        'binary': {
            'hipp_z2.0': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='Hippocampus (z<-2.0)'),
            'hipp_z2.5': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='Hippocampus (z<-2.5)'),
            'mttvol_z2.0': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False, nickname='MTV (z<-2.0)'),
            'mttvol_z2.5': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False, nickname='MTV (z<-2.5)')},
        'categorical': {
            'hipp_quantiles': Quantiles('HIPPOCAMPUS_VOL', nickname='Hippocampus (Quantiles)'),
            'mttvol_quantiles': Quantiles('META_TEMPORAL_VOL', nickname='MTV (Quantiles)')},
        'continuous': {
            'hipp': Continuous('HIPPOCAMPUS_VOL', nickname='Hippocampus'),
            'mttvol': Continuous('META_TEMPORAL_VOL', nickname='MTV')}},
    }

BINARY_DATA_DRIVEN = {
    'amyloid': {
        'binary': {
            'composite_gmm': BinaryGMM('AMYLOID_COMPOSITE'),
            'centiloid_gmm': BinaryGMM('Centiloid'),
            'composite_z2.0': BinaryZScore('AMYLOID_COMPOSITE', zcutoff=2.0, control_col='Control'),
            'composite_z2.5': BinaryZScore('AMYLOID_COMPOSITE', zcutoff=2.5, control_col='Control'),
            'centiloid_z2.0': BinaryZScore('Centiloid', zcutoff=2.0, control_col='Control'),
            'centiloid_z2.5': BinaryZScore('Centiloid', zcutoff=2.5, control_col='Control')}
        },
    'tau': {
        'binary': {
            'mtt_gmm': BinaryGMM('META_TEMPORAL_TAU'),
            'mtt_z2.0': BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.0),
            'mtt_z2.5': BinaryZScore('META_TEMPORAL_TAU', 'Control', zcutoff=2.5)}
        },
    'neurodegeneration': {
        'binary': {
            'hipp_z2.0': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False),
            'hipp_z2.5': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False),
            'mttvol_z2.0': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False),
            'mttvol_z2.5': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False)},
        }
    }

BINARY_ESTABLISHED = {
    'amyloid': {
        'binary': {
            'composite_1.11': BinaryManual('AMYLOID_COMPOSITE', 1.11),
            'centiloid_20': BinaryManual('Centiloid', 20)}
        },
    'tau': {
        'binary': {
            'mtt_1.20': BinaryManual('META_TEMPORAL_TAU', cutoff=1.20),
            'mtt_1.21': BinaryManual('META_TEMPORAL_TAU', cutoff=1.21),
            'mtt_1.23': BinaryManual('META_TEMPORAL_TAU', cutoff=1.23),
            'mtt_1.33':  BinaryManual('META_TEMPORAL_TAU', cutoff=1.33)},
        },
    'neurodegeneration': {
        'binary': {
            'hipp_z2.0': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.0, greater=False),
            'hipp_z2.5': BinaryZScore('HIPPOCAMPUS_VOL', control_col='Control', zcutoff=-2.5, greater=False),
            'mttvol_z2.0': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.0, greater=False),
            'mttvol_z2.5': BinaryZScore('META_TEMPORAL_VOL', control_col='Control', zcutoff=-2.5, greater=False)},
        },
    }
