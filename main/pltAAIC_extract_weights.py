
# ------ imports ------

import os

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from common import load_results

# ------ load models ------

models = load_results('expaaic_svms', 'models.pickle')

# ------ define helpers ------

def chen_2023_haufe_weights(train_X, y_pred):
    # Code is taken from JianZhong Chen's 2023 NeuroImage paper
    # https://doi.org/10.1016/j.neuroimage.2023.120115
    # https://github.com/ThomasYeoLab/CBIG/blob/424e2d9633521b2521d46e4bf22efa74d894f29f/stable_projects/predict_phenotypes/ChenOoi2023_ICCW/regression/RF/CBIG_ICCW_RF_interpretation.py#L276-L283
    # There, it was used to estimate Haufe feature importances
    # for nonlinear RF models.

    # The code is used mostly as is, except we don't redo
    # demeaning of X_train (this is already done for the model training)

    # Took me some time to understand this code, but think I get it now.
    # Covariance formula is as follows:
    #   cov(X, Y) = E[(X-E[X]) * (Y-E[Y])]

    # X-E[X] is the variable demeaned, which is achieved in the first two
    # lines of code for the training data and targets

    # The outer expectation is then the *mean* of the product of the
    # demeaned X and Y.  The last two lines achieve this, by
    # taking the dot product of each feature with Y and then dividing
    # by the number of observations (train_X.shape[0]).  The dot product
    # sums all the products, and the division takes the average.

    # Basically, we get the covariance between each feature and the model
    # predictions.  In Haufe 2014 (https://doi.org/10.1016/j.neuroimage.2013.10.067),
    # this is described in the "Simplifying conditions" section.  We can
    # take advantage of this since we only have one target (K=1).

    demean_y_pred = y_pred - np.mean(y_pred)
    demean_X_train = train_X.T

    fivals_Haufe = np.dot(demean_X_train, demean_y_pred) / train_X.shape[0]
    fivals_Haufe = fivals_Haufe[:, np.newaxis]

    return fivals_Haufe

def get_haufe_weights(dataset, models, repeats, outer_splits, outer_seed, stratify):
    feats = len(models[0].predictors)
    Amat = np.zeros((feats, len(models)))
    # Wmat = np.zeros((feats, len(models)))
    for r in range(repeats):
        outer_cv = StratifiedKFold(n_splits=outer_splits, random_state=outer_seed + r, shuffle=True)
        # outer CV loop
        for i, (outer_train_index, outer_test_index) in enumerate(outer_cv.split(dataset, dataset[stratify])):
            outer_train = dataset.iloc[outer_train_index, :]
            idx = (r * outer_splits) + i
            model = models[idx]

            # get weight-independent Haufe importances
            # X_unscaled is needed for predictions
            # X_scaled is needed for feature importances
            X_unscaled = outer_train[model.predictors].values
            X_scaled = model.pipeline['scaler'].transform(X_unscaled)
            Y_pred = model.pipeline.predict(X_unscaled)
            A = chen_2023_haufe_weights(X_scaled, Y_pred)
            Amat[:, idx] = A[:, 0]

            # previous approach to get importances from linear SVM
            # only possible with linear kernel
            # as you can get a coefficient for each feature
            # X_unscaled = outer_train[model.predictors].values
            # X_scaled = model.pipeline['scaler'].transform(X_unscaled)
            # W = model.pipeline['svm'].coef_[np.newaxis, :]
            # Y = outer_train[model.target].values[:, np.newaxis]
            # A = haufe_trick(W, X_scaled, Y)
            # Amat[:, idx] = A[:, 0]
            # Wmat[:, idx] = W

    Adf = pd.DataFrame(Amat.T, columns=model.predictors)
    return Adf

# ------ define helpers ------
dataset = pd.read_csv('datasets/maindata.csv')

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

outputfolder = os.path.join('outputs', 'haufe_weights_aaic')
if not os.path.exists(outputfolder):
    os.mkdir(outputfolder)

for name, model_list in models.items():
    print(f'{name}...')
    fixname = ''.join([c for c in name if c not in [' ', '[', ']']])
    haufe_out = os.path.join(outputfolder, f'svm_weights_haufe_{fixname}.csv')
    A = get_haufe_weights(dataset=dataset,
                          models=model_list,
                          repeats=repeats,
                          outer_splits=outer_splits,
                          outer_seed=outer_seed,
                          stratify=stratify)
    A.to_csv(haufe_out, index=False)

print('Done!')

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
