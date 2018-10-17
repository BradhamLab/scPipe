import pandas as pd
import numpy as np
import scanpy.api as sc

from sklearn import mixture
from scipy import linalg

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

def exp_gaus_mixture(X):
    """
    Fit data to a mixture of an exponential and gaussian distribution.
    
    Parameters
    ----------
    X : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """

    #  use count data for mixtures
    counts, bins, patches = plt.hist(X, bins=30)
    mode_idx = np.where(counts == np.max(counts[1:]))[0][0]
    normal_max = bins[-1]
    normal_min = bins[mode_idx - len(bins[mode_idx:])]
    sd = np.std(X[X >= normal_min])
    mean = np.mean(X[X >= normal_min])

    return 0