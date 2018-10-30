import pandas as pd
import numpy as np
import scanpy.api as sc

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

    # read in data expression data
    expr_data = pd.read_csv(expr_file, index_col=0)
    if transpose_expr:
        expr_data = expr_data.T

    # read in gene annotations, fix ids between annotations + alignments
    gene_data = pd.read_csv(gene_file, index_col=0)
    gene_data['scaffold'] = [x.replace('model', 'TU') for x in gene_data.index]
    gene_data = gene_data.set_index('scaffold')
    
    # read in cell annotations
    cell_data = pd.read_csv(cell_file, index_col=0)

    # add null cell annotations for missing cells
    null_cells = list(set(expr_data.index).difference(cell_data.index))
    if len(null_cells) > 0:
        null_cell_anno = pd.DataFrame(index=null_cells,
                                      columns=cell_data.columns)
        cell_data = pd.concat([cell_data, null_cell_anno])

    # remove specified cells for filtering
    keep_cells = cell_data.index.values
    if filter_cells is not None:
        keep_cells = [x for x in keep_cells if x not in filter_cells]
        
    cell_data = cell_data.filter(keep_cells, axis=0)
    expr_data = expr_data.filter(keep_cells, axis=0)

    # set null annotations to unannotated genes
    null_genes = list(set(expr_data.columns).difference(gene_data.index))
    if len(null_genes) > 0:
        null_gene_anno = pd.DataFrame(index=null_genes,
                                      columns=gene_data.columns)
        gene_data = pd.concat([gene_data, null_gene_anno])

    # filter genes to only include genes with reads
    gene_data = gene_data.filter(expr_data.columns, axis=0)

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


def gini_coefficient(X):
    """
    Calculate the Gini coefficient of a given expression profile.

    Calculates the Gini coefficient to measure dispersion of counts. Calculated
    using the median absolute deviation form.
    
    Parameters
    ----------
    X : np.array
        Expression profile for a given gene. 
    
    Returns
    -------
    float
        The Gini coefficient to measure "inequality" within the profile.
    
    References
    ----------
         Damgaard, Christian. "Gini Coefficient." From MathWorld--A Wolfram Web
         Resource, created by Eric W. Weisstein.
         http://mathworld.wolfram.com/GiniCoefficient.html 
    """

    return np.sum(np.abs(np.subtract.outer(X, X))) / (2*len(X)*np.sum(X))


def variable_genes(anno_df, method='gini', percentile=0.75, ignore_zeros=True):
    """
    Return most variable genes, measured by dispersion. 
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of gene expression data. 
    method : str, optional
        Method to calculate measure dispersion/variance within an expression
        profile. Default is 'gini', and the Gini coefficient will be calculated.
        Possible values include: 'gini', 'dispersion'. 
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
        List of ordered gene ids from least to most variable. 
    """
    gene_dispersion = np.zeros(anno_df.shape[1])
    if not 0 <= percentile <= 1:
        raise ValueError('`percentile` must be a value in between 0 and 1.')
    if not ignore_zeros:
        Warning("Keeping zeros when calculating dispersion often leads to NaNs.")
    if method not in ['gini', 'dispersion']:
        raise ValueError('Unsupported dispersion method: {}'.format(method))
    for i, gene in enumerate(anno_df.var.index):
        expression = np.array(anno_df[:, gene].X)
        if ignore_zeros:
            expression = expression[expression != 0]
        if len(expression) > 0:
            if method == 'gini':
                gene_dispersion[i] = gini_coefficient(expression)
            elif method == 'dispersion':
                gene_dispersion[i] = dispersion(expression)
    sorted_dispersion = anno_df.var.index.values[np.argsort(gene_dispersion)]
    start_idx = min(int(len(sorted_dispersion)*percentile),
                    anno_df.shape[1] - 1)

    return sorted_dispersion[start_idx:]


def filter_by_coverage(anno_df, threshold=0.5):
    """
    Filter genes below a specified dropout percent. 
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of expression data.
    threshold : float, optional
        Minimum percent of non-zeros to needed to keep an expression profile
        from a certain gene. The default is 0.5, which will filter genes with
        zeros for greater than 50% of measurements.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with low-coverage genes removed.
    """
    coverage = np.array([sum(anno_df.X[:, i] != 0) / anno_df.shape[0] for i in \
                         range(anno_df.shape[1])])
    keep_genes = anno_df.var.index.values[np.where(coverage >= threshold)]
    return anno_df[:, keep_genes]


def correlated_genes(anno_df, genes, threshold, mask_zeros=False,
                     min_pairwise=None):
    """
    Find genes correlated with a set of provided genes.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe with expression profiles. 
    genes : list, str
        List of gene indices to check for correlation relationships.
    threshold : float
        $R^2$ threshold for correlation considerations.
    mask_zeros : [type], optional
        Whether to consider zeros when calculating correlations.
    min_pairwise : int, optional
        Minimum number of observations with pairwise values required for a
        correlation to be considered valid. Used when `mask_zeros=True`. Default
        in None.
    
    Returns
    -------
    pd.DataFrame
        Dataframe with correlation results and gene annotations. 
    """

    expr = pd.DataFrame(data=anno_df.X, columns=anno_df.var.index.values,
                        index=anno_df.obs.index.values)
    if mask_zeros:
        expr[expr == 0] = None

    corr_mat = expr.corr(min_periods=min_pairwise)
    result_dict = {'Subset.Gene': [],
                   'Correlated.Gene': [],
                   'R2': [],
                   'R': [],
                   'Percent.Zero': []}
    all_genes = set(genes)
    for x in genes:
        for y in corr_mat.columns.values:
            r_square = corr_mat.loc[x, y]**2
            if r_square >= threshold and r_square != 1:
                result_dict['Subset.Gene'].append(x)
                result_dict['Correlated.Gene'].append(y)
                result_dict['R2'].append(r_square)
                result_dict['R'].append(corr_mat.loc[x, y])
                result_dict['Percent.Zero'].append(sum(anno_df[:, y].X == 0) / \
                                                    anno_df.shape[0])
                all_genes.add(y)
    if mask_zeros:
        anno_df.X[np.isnan(anno_df.X)] = 0

    results = pd.DataFrame(result_dict)
    results['Subset.Name'] = results.apply(lambda x:
                              anno_df.var.loc[x['Subset.Gene'], 'UniProt.Name'],
                              axis=1)
    results = pd.merge(results, anno_df.var, left_on=results['Correlated.Gene'],
                   right_on=anno_df.var.index.values)
    results = results.drop('key_0', axis=1)
    return results