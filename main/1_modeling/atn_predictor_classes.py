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
    
    def __init__(self, atn=None, variable_type=None, nickname=None):
        self.atn = atn
        self.variable_type = variable_type
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

    def __init__(self, y_col, control_col, zcutoff=2.0, greater=True, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='binary')
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

    def __init__(self, y_col, greater=True, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='binary')
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

    def __init__(self, y_col, cutoff, greater=True, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='binary')
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
        regions_in_stage = (np.array(groupings) == i)
        n_regions_in_stage = regions_in_stage.sum()
        sub = data[:, regions_in_stage]
        freqs = sub.sum(axis=1) / n_regions_in_stage

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

    def __init__(self, y_col, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='categorical')
        self.y_col = y_col

        self.quantiles = None

    def fit(self, data):
        self.quantiles = np.quantile(data[self.y_col], [.25, .5, .75, 1])

    def covariates(self, data):
        digitized = np.digitize(data[self.y_col], self.quantiles)
        return pd.get_dummies(digitized).to_numpy().astype(float)

class CategoricalStager(ATNPredictor):

    def __init__(self, columns, groupings=None, method='gmm', non_stageable='NS',
                 p='any', atn=None, nickname=None, **kwargs):
        super().__init__(atn=atn, nickname=nickname, variable_type='categorical')
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
        return pd.get_dummies(stages).to_numpy().astype(float)
    
class GMMWithIndeterminateZone(ATNPredictor):

    def __init__(self, y_col, margin=.1, greater=True, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='categorical')
        self.y_col = y_col
        self.margin = margin
        self.greater = greater
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
        if self.margin == 0:
            return np.where(self.operator(y, self.cutoff), 1., 0.)[:, np.newaxis]
        
        dim = 1 if self.greater else 0
        probs = self.gmm.predict_proba(y.values.reshape(-1, 1))[:, dim]
        digitized = np.digitize(probs, [0.5 - self.margin, 0.5 + self.margin])
        return pd.get_dummies(digitized).to_numpy().astype(float)

# Continuous predictors ----

class Continuous(ATNPredictor):

    def __init__(self, y_col, atn=None, nickname=None):
        super().__init__(atn=atn, nickname=nickname, variable_type='continuous')
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
            if 'gamma' in kwargs:
                del kwargs['gamma']
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
