
# ------ imports ------

import os
import pickle
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.covariance import EmpiricalCovariance

# ------ load models ------
# make sure ATN predictor classes are avaiable
sys.path.append('../1_modeling')

this_dir = os.path.dirname(os.path.abspath(__file__))
exp1_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp1_svms_global_cognition'))
exp1_results = os.path.join(exp1_folder, 'models.pickle')
    
if not os.path.exists(exp1_results):
    short_folder = os.path.abspath(os.path.join(this_dir, '..', '..', 'outputs', 'exp1_svms_global_cognition_short'))
    short_results = exp1_results = os.path.join(short_folder, 'models.pickle')
    
    if os.path.exists(short_results):
        warnings.warn('The pickled models are missing for exp1 (SVMs), but the short results are present.  Using the latter!',
                      RuntimeWarning)
        exp1_folder = short_folder
        exp1_results = short_results
    else:
        raise RuntimeError('Results for exp1 are missing.')

with open(exp1_results, 'rb') as f:
    models = pickle.load(f)

# ------ define helpers ------

# https://gist.github.com/wmvanvliet/08f9b4b8c5e4f1bb0a25d5eceae1b419
def haufe_trick(W, X, Y):
    """Perform the Haufe trick."""
    # Computing the covariance of X and Y involves removing the mean
    cov_X = EmpiricalCovariance().fit(X).covariance_
    cov_Y = EmpiricalCovariance().fit(Y).covariance_

    # The Haufe trick
    A_hat = cov_X.dot(W.T).dot(np.linalg.pinv(cov_Y))
    return A_hat

def get_haufe_weights(dataset, models, repeats, outer_splits, outer_seed, stratify):
    feats = len(models[0].predictors)
    Amat = np.zeros((feats, len(models)))
    Wmat = np.zeros((feats, len(models)))
    for r in range(repeats):
        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)
        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(dataset, dataset[stratify])):
            outer_train = dataset.iloc[outer_train_index, :]
            idx = (r * outer_splits) + i
            model = models[idx]
            X_unscaled = outer_train[model.predictors].values
            X_scaled = model.pipeline['scaler'].transform(X_unscaled)
            W = model.pipeline['svm'].coef_[np.newaxis, :]
            Y = outer_train[model.target].values[:, np.newaxis]
            A = haufe_trick(W, X_scaled, Y)
            Amat[:, idx] = A[:, 0]
            Wmat[:, idx] = W
    
    Wdf = pd.DataFrame(Wmat.T, columns=model.predictors)
    Adf = pd.DataFrame(Amat.T, columns=model.predictors)
    return Wdf, Adf
    
# ------ define helpers ------
dataset = pd.read_csv('../../outputs/maindata/maindata.csv')

# ------ SVM parameters ------
# recreate sampling loop from SVM experiment
# IMPORTANT !!!
# These parameters should match the SVM experiment file
# otherwise, different data splitting will occur
repeats=10
outer_splits=10
outer_seed=0
stratify='CDRBinned'

# ------ Run! ------
        
for name, model_list in models.items():
    fixname = ''.join([c for c in name if c not in [' ', '[', ']']])
    raw_out = os.path.join(exp1_folder, f'svm_weights_raw_{fixname}.csv')
    haufe_out = os.path.join(exp1_folder, f'svm_weights_haufe_{fixname}.csv')
    W, A = get_haufe_weights(dataset=dataset,
                             models=model_list,
                             repeats=repeats,
                             outer_splits=outer_splits,
                             outer_seed=outer_seed,
                             stratify=stratify)
    W.to_csv(raw_out, index=False)
    A.to_csv(haufe_out, index=False)

# ------ sanity check plots for Haufe transform ------
# import matplotlib.pyplot as plt
# from matplotlib import colors

# # svm weights
# data = Wbig
# mag = np.abs(data).max()
# plt.figure(figsize=(10, 10))
# plt.imshow(data, cmap='RdBu',
#            norm=colors.TwoSlopeNorm(vmin=-mag,
#                                     vcenter=0.,
#                                     vmax=+mag))
# plt.xlabel('Iteration')
# plt.ylabel('Feature')
# plt.title('Raw SVM weights')
# plt.colorbar(fraction=0.046, pad=0.04)
# plt.tight_layout()
# plt.savefig('/home/tom.earnest/Desktop/svm_weights.png')

# # haufe weights
# data = Abig
# mag = np.abs(data).max()
# plt.figure(figsize=(10, 10))
# plt.imshow(data, cmap='RdBu',
#            norm=colors.TwoSlopeNorm(vmin=-mag,
#                                     vcenter=0.,
#                                     vmax=+mag))
# plt.xlabel('Iteration')
# plt.ylabel('Feature')
# plt.title('Haufe-transformed SVM weights')
# plt.colorbar(fraction=0.046, pad=0.04)
# plt.tight_layout()
# plt.savefig('/home/tom.earnest/Desktop/haufe_weights.png')

# # cov mat
# X_unscaled = outer_train[model.predictors].values
# X_scaled = model.pipeline['scaler'].transform(X_unscaled)
# data = EmpiricalCovariance().fit(X_scaled).covariance_
# mag = np.abs(data).max()
# plt.figure(figsize=(10, 10))
# plt.imshow(data, cmap='RdBu',
#            norm=colors.TwoSlopeNorm(vmin=-mag,
#                                     vcenter=0.,
#                                     vmax=+mag))
# plt.xlabel('Feature')
# plt.ylabel('Feature')
# plt.title('Covariance matrix')
# plt.colorbar(fraction=0.046, pad=0.04)
# plt.tight_layout()
# plt.savefig('/home/tom.earnest/Desktop/covariance.png')


