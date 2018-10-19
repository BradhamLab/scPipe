import pandas as pd
import numpy as np
import scanpy.api as sc

from skimage.filters import threshold_otsu
from scipy import linalg, stats
import pomegranate as pmg

from matplotlib import pyplot as plt

def get_gene_identifier(scaffold, gene_df):
    """
    Get gene name associated with scaffold id.
    
    Parameters
    ----------
    scaffold : str
        Scaffold transcript id in `gene_df` as an index value.
    gene_df : pd.DataFrame
        A (gene x feature) dataframe with 'UniProt.Name', 'UniProt.ID', 'SPU',
        and 'NCBI.ID' columns containing different gene id/name values.
    
    Returns
    -------
    str
        Gene name associated with `scaffold` id. Name priority follows
        'UniProt.Name', 'Uniprot.ID', 'SPU', 'NCBI.ID', `scaffold` in order. 
    """

    if not pd.isnull(gene_df.loc[scaffold, 'UniProt.Name']):
        return gene_df.loc[scaffold, 'UniProt.Name'].split('_')[0]
    elif not pd.isnull(gene_df.loc[scaffold, 'UniProt.ID']):
        return gene_df.loc[scaffold, 'UniProt.ID']
    elif not pd.isnull(gene_df.loc[scaffold, 'SPU']):
        return gene_df.loc[scaffold, 'SPU']
    elif not pd.isnull(gene_df.loc[scaffold, 'NCBI.ID']):
        return gene_df.loc[scaffold, 'NCBI.ID']
    else:
        return scaffold


def create_annotated_df(expr_file, gene_file, cell_file, filter_cells=None,
                        transpose_expr=True):
    """
    Create an AnnData object to hold expression, cell, and gene data.
    
    Parameters
    ----------
    expr_file : string
        File path to csv file of expression data. Should be gene x cell.
    gene_file : string
        File path to gene annotation file. Should be gene x feature.
    cell_file : string
        File path to cell annotation file. Should be cell x feature. 
    filter_cells : list, optional
        Cells to remove from the final AnnData (the default is None, which
        will not filter cells.)
    tranpose_expr : boolean, optional
        Whether to transpose the expression matrix. Set true if file is written
        as (genes x sample) instead of (sample x gene). Default is True.
    
    Returns
    -------
    sc.AnnData
        Annotated data frame containing expression data, cell data, and gene
        data.
    """

    # read in data
    expr_data = pd.read_csv(expr_file, index_col=0)
    if transpose_expr:
        expr_data = expr_data.T

    gene_data = pd.read_csv(gene_file, index_col=0)
    gene_data['scaffold'] = [x.replace('model', 'TU') for x in gene_data.index]
    gene_data = gene_data.set_index('scaffold')  # this is filling with NaN

    cell_data = pd.read_csv(cell_file, index_col=0)

    # only keep cells found in both cell data and expr data
    keep_cells = list(set(cell_data.index).intersection(expr_data.index))

    # remove non-PMC cells
    if filter_cells is not None:
        keep_cells = [x for x in keep_cells if x not in filter_cells]
        
    cell_data = cell_data.filter(keep_cells, axis=0)
    expr_data = expr_data.filter(keep_cells, axis=0)

    # remove genes that are not in both annotation and and expression datasets 
    keep_genes = list(set(expr_data.columns).intersection(gene_data.index))
    gene_data = gene_data.filter(keep_genes, axis=0)
    expr_data = expr_data.filter(keep_genes, axis=1)

    # create anno df formatted cells x genes 
    anno_df = sc.AnnData(expr_data.values, obs=cell_data, var=gene_data)
    return anno_df


def filter_cells(anno_data, cell_col, filter_func):
    """
    Filter an annotated data frame to  cells only.
    
    Parameters
    ----------
    anno_data : sc.AnnData
        Annotated data frame with cell annotations denoting treatment.
    cell_col : string
        Name of column in `sc.AnnData.obs` containing features to filter on.
    filter_func : function
        Function that takes values from `cell_col`, and returns a boolean value
        whether a condition is met. Rows where the function returns True will
        be kept, while rows evaluated as False will be removed. 
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with filtered cells.
    """
    keep_cells = [x for x in anno_data.obs.index if\
                  filter_func(anno_data.obs.loc[x, cell_col])]
    return anno_data[keep_cells, ]


def filter_genes(anno_data, gene_col, filter_func):
    """
    Filter an annotated data frame to  cells only.
    
    Parameters
    ----------
    anno_data : sc.AnnData
        Annotated data frame with cell annotations denoting treatment.
    gene_col : string
        Name of column in `sc.AnnData.var` containing features to filter on.
    filter_func : function
        Function that takes values from `gene_col`, and returns a boolean value
        whether a condition is met. Rows where the function returns True will
        be kept, while rows evaluated as False will be removed.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with filtered genes.
    """
    keep_genes = [x for x in anno_data.var.index if\
                  filter_func(anno_data.var.loc[x, gene_col])]
    return anno_data[:, keep_genes]


def dispersion(array):
    """
    Measure dispersion within a numerical array.
    
    Parameters
    ----------
    array : np.ndarray
        Array of gene expression values.
    
    Returns
    -------
    float
        Gene expression dispersion measured by (Q3 - Q1) / (Q3 + Q1), where Q3
        is the third quartile and Q1 is the first quartile.

    References
    ----------
        Wiki:
        https://en.wikipedia.org/wiki/Quartile_coefficient_of_dispersion

        Original Paper:
        Bonett, D. G. (2006). "Confidence interval for a coefficient of quartile
        variation". Computational Statistics & Data Analysis. 50 (11): 
        2953–2957. 
    """

    quartiles = np.percentile(array, [25, 75])
    return (quartiles[1] - quartiles[0]) / (quartiles[1] + quartiles[0])


def variable_genes(anno_df, percentile=0.75, ignore_zeros=True):
    """
    Return most variable genes, measured by dispersion. 
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of gene expression data. 
    percentile : float, optional
        Percentile of most variable genes to return. Value should be between
        0 and 1. Default is 0.75, and only the top 25% of genes will be
        returned. 
    ignore_zeros : bool, optional
        Whether to ignore zeros during dispersion calculation. Default is True.
    
    Raises
    ------
    ValueError
        Raised if `percentile` is not between 0 and 1. 
    
    Returns
    -------
    numpy.ndarray
        List of gene ids of most variable genes. 
    """
    gene_dispersion = np.zeros(anno_df.shape[1])
    if not 0 < percentile < 1:
        raise ValueError('`percentile` must be a value in between 0 and 1.')
    if not ignore_zeros:
        Warning("Keeping zeros when calculating dispersion often leads to NaNs.")
    for i, gene in enumerate(anno_df.var.index):
        expression = anno_df[:, gene].X
        if ignore_zeros:
            expression = expression[expression != 0]
        if len(expression) > 0:
            gene_dispersion[i] = dispersion(expression)
    sorted_dispersion = np.argsort(gene_dispersion)
    start_idx = int(len(sorted_dispersion)*percentile)

    return anno_df.var.index.values[start_idx:]


def binarize_expression(anno_df, filter_zeros):
    """
    Binarize gene expression in cells by fitting gaussian mixture models.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of expression data.
    filter_zeros : boolean
        Whether to remove zeros before fitting mixture model.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with expression data set to 1 or 0 depending whether
        the gene is 'on' or 'off', respectively, in each cell.
    """



    return None

def noise_signal_mixture(X):
    """
    Fit data to a mixture model to identify signal and technical noise.

    Fit data to a mixture model where technical noise of scRNAseq profiles is
    modelled using a Poisson distribution, and true expression is modelled as a
    Gaussian Distribution. 
    
    Parameters
    ----------
    X : numpy.array
        single-cell expression profile across cells.
    
    Returns
    -------
    pomegranate.GeneralMixtureModel
        General mixture model with two components: a Poisson distribution to
        model technical noise, and a Gaussian distribution to model signal.
    """

    #  use count data for mixtures
    counts, bins = np.histogram(X, bins=30)
    mode_idx = np.where(counts == np.max(counts[1:]))[0][0]
    normal_spread = bins[-1] - bins[mode_idx]
    normal_min = bins[np.where(bins < bins[mode_idx] - normal_spread)[0][-1]]

    signal = pmg.NormalDistribution.from_samples(X[X >= normal_min])

    pois_lambda = max(np.mean(X[X < normal_min]), sum(X==0)/len(X))

    noise = pmg.PoissonDistribution(pois_lambda)
    gmm = pmg.GeneralMixtureModel([noise, signal])
    gmm = gmm.fit(X)

    return gmm


def threshold_expression(X, method='otsu'):
    """
    Threshold gene expression using pre-defined methods.
    
    Parameters
    ----------
        X: np.array
            Expression profile across cells.
    
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

    threshold = 0
    if method == 'otsu':
        threshold = threshold_otsu(X)
    elif method == 'mixture':
        gmm = noise_signal_mixture(X)
        space = np.arange(0, np.max(X), 0.1).reshape(-1, 1)
        probabilities = gmm.predict_proba(space)
        idx = 0
        while probabilities[idx][0] > probabilities[idx][1] and \
        idx < probabilities.shape[0]:
            idx += 1
        threshold = space[idx][0]
    else:
        raise ValueError("Unsupported method: {}".format(method))
    X[X < threshold] = 0
    return (X, threshold)
