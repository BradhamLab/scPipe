import numpy as np
from scipy import stats
from anndata.base import ArrayView
from skimage.filters import threshold_otsu
import re

import sc_classes
import sc_utils

def set_on_off(anno_df, method='mixture', test_fits=True, overwrite=False):
    """
    Set genes to an on or off state by fitting a general mixture model.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of expression data.
    method : str, optional
        Method to determine whether a gene is 'on' or 'off' in a given cell.
        Default is 'mixture', which will fit a Poisson-Gaussian mixture model to
        distinguish between technical noise and actual signal. See 
        `threshold_expression()` for all methods.
    test_fits : boolean, optional
        If `method == True`, whether to test mixture model against an all
        signal model (Gaussian) and an all noise model (Poisson). Fits are
        compared using Kamalgorov-Smirnoff scores. If an all noise model
        performs better, all counts will be set to zero. If the all signal model
        performs better, all counts are kept. If the mixture fits the
        emperical distribution best, Poisson sample counts are set to zero.
        Default is True.
    overwrite : boolean, optional
        Wether to write binarized expression values to original `sc.AnnData`
        object. Default is `False`, where a copy of the object will be created
        and original counts will be unmolested. 
    
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with expression data set to 1 or 0 depending whether
        the gene is 'on' or 'off', respectively, in each cell.
    """
    binarized = anno_df
    if not overwrite:
        binarized = anno_df.copy()
    for gene in binarized.var.index:
        if method == 'mixture':
            if test_fits:
                out = fit_expression(binarized[:, gene].X)
            else:
                out, __ = threshold_expression(binarized[:, gene].X,
                                                method='mixture')
        elif method == 'otsu':
            out, __ = threshold_expression(binarized[:, gene].X, method='otsu')
        binarized[:, gene].X = out
    return binarized
        

def fit_expression(X, silent=True):
    mu = np.mean(X)
    sigma = np.std(X)

    noise = stats.poisson(mu)
    signal = stats.norm(loc=mu, scale=sigma)
    mixture = sc_classes.SignalNoiseModel(X)
    models = [noise.cdf, signal.cdf, mixture.cdf]
    scores = [stats.kstest(X, f)[0] for f in models]

    out = X.copy()
    if scores[0] == np.min(scores):
        msg = 'All reads likely noise. Setting counts to zero.'
        out = np.zeros_like(X)
    elif scores[1] == np.min(scores):
        msg = 'All reads likely signal. Keeping counts.'
    else:
        msg = 'Reads likely a mixture of noise and signal. Setting noise reads'\
        +     ' to zeros.'
        out, __ = threshold_expression(out, value=mixture.threshold())
    if not silent:
        print(msg)
    return out


def threshold_expression(X, value=None, method='otsu'):
    """
    Threshold gene expression using pre-defined methods.
    
    Parameters
    ----------
        X: np.array
            Expression profile across cells.
        
        value: float, optional
            Value to threshold expression by. Default is `None`, and a threshold
            will be estimated using the method passed in the `method` argument.
    
        method: str, optional
             Method to determine expression threshold. Default is otsu, which
             will perform Otsu thresholding. Additional method is 'mixture',
             which fits a mixture model to differentiate technical noise from
             signal.

    Returns
    -------
        tuple, (np.array, float)
            Tuple of thresholded expression profile with filtered values set to
            zero and discovered threshold.
    """
    if isinstance(X, ArrayView):
        X = np.array(X)
    threshold = 0
    if value is not None:
        threshold = value
    elif method == 'otsu':
        threshold = threshold_otsu(X)
    elif method == 'mixture':
        gmm = sc_classes.SignalNoiseModel(X)
        threshold = gmm.threshold()
    else:
        raise ValueError("Unsupported method: {}".format(method))
    X[X < threshold] = 0
    return (X, threshold)


def compare_non_zero_expression(anno_df, gene, group_col, groups,
                               compare_zeros=True):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    gene : [type]
        [description]
    group_col : [type]
        [description]
    groups : [type]
        [description]
    compare_zeros : bool, optional
        [description] (the default is True, which [default_description])
    
    Returns
    -------
    [type]
        [description]
    """

    regex = [re.compile(x) for x in groups]
    group_1 = sc_utils.filter_cells(anno_df, group_col,
                                    lambda x: regex[0].search(x) is not None)

    group_2 = sc_utils.filter_cells(anno_df, group_col,
                                    lambda x: regex[1].search(x) is not None)
    # create output dictionary
    out_dict = {'group1': groups[0], 'group2': groups[1]}
    # just bake other fisher stuff into this
    non_zero_g1 = group_1[:, gene].X[group_1[:, gene].X != 0]
    non_zero_g2 = group_2[:, gene].X[group_2[:, gene].X != 0]

    if compare_zeros:
        g1_not = len(non_zero_g1)
        g1_zeros = group_2.shape[0] - g1_not
        g2_not = len(non_zero_g2)
        g2_zeros = group_2.shape[0] - g2_not
        #                   group1    group 2
        counts = np.array([[g1_zeros, g2_zeros], # zero expression
                           [g1_not,   g2_not]])  # non-zero expression
        odds, fisher_p = stats.fisher_exact(counts)
        out_dict['fisher.odds'] = odds
        out_dict['fisher.pvalue'] = fisher_p
        out_dict['group1.countZero'] = counts[0, 0]
        out_dict['group2.countZero'] = counts[0, 1]
        out_dict['group1.countNot'] = counts[1, 0]
        out_dict['group2.countNot'] = counts[1, 1]
    t_val, t_p = stats.ttest_ind(a=non_zero_g1, b=non_zero_g2)
    out_dict['t.value'] = t_val
    out_dict['ttest.pvalue'] = t_p
    out_dict['group1.avg'] = np.mean(non_zero_g1)
    out_dict['group2.avg'] = np.mean(non_zero_g2)
    out_dict['group1.std'] = np.std(non_zero_g1)
    out_dict['group2.std'] = np.std(non_zero_g2)
    out_dict['gene'] = gene

    return out_dict
