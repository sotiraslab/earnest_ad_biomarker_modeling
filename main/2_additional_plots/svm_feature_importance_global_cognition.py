
import pickle
import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

# load models
# make sure ATN predictor classes are avaiable
sys.path.append('../1_modeling')
models_path = '../../outputs/exp1_svms_global_cognition_short/models.pickle'
with open(models_path, 'rb') as f:
    models = pickle.load(f)
    
amy_svm = models['Amyloid SVM']
tau_svm = models['Tau SVM']
gm_svm = models['GM SVM']

# define haufe transform

# https://gist.github.com/wmvanvliet/08f9b4b8c5e4f1bb0a25d5eceae1b419
def haufe_trick(W, X, Y):
    """Perform the Haufe trick."""
    # Computing the covariance of X and Y involves removing the mean
    X_ = X - X.mean(axis=0)
    Y_ = Y - Y.mean(axis=0)
    cov_X = X_.T.dot(X_)
    cov_Y = Y_.T.dot(Y_)

    # The Haufe trick
    A_hat = cov_X.dot(W.T).dot(np.linalg.pinv(cov_Y))
    return A_hat

# load data
dataset = pd.read_csv('../../outputs/maindata/maindata.csv')

# recreate sampling loop from SVM experiment
# IMPORTANT !!!
# These parameters should match the SVM experiment file
# otherwise, different data splitting will occur
repeats=10
outer_splits=10
outer_seed=0
stratify='CDRBinned'

# loop
Abig = np.zeros((85, 100))
for r in range(repeats):
    outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)
    # outer CV loop
    for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(dataset, dataset[stratify])):
        outer_train = dataset.iloc[outer_train_index, :]
        
        idx = (r * outer_splits) + i
        model = tau_svm[idx]
        X = outer_train[model.predictors].values
        W = model.pipeline['svm'].coef_[np.newaxis, :]
        Y = outer_train[model.target].values[:, np.newaxis]
        A = haufe_trick(W, X, Y)
        Abig[:, idx] = A[:, 0]
        


