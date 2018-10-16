import string

import igraph
import louvain
import matplotlib.pyplot as plt
from matplotlib.text import OffsetFrom
import numpy as np
import pandas as pd
import scanpy.api as sc
import seaborn as sns
import sklearn
import umap
from scipy import spatial, stats

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=80)
sc.logging.print_versions()

expr_file = '../../output/final/normalized_tpm_matrix.csv'
gene_file = '/home/dakota/SequenceData/evm_annotations.csv'
cell_file = '../../output/metadata/filtered_metadata.csv'

bad_cells = ["ASW_C07_2016-07-29", "ASW_C01_2016-07-29",
             "ASW_C09_2016-07-29", "ASW_A02_2016-07-29",
             "ASW_G06_2016-07-29", "ASW_E05_2016-07-29",
             "ASW_A12_2016-07-29", "ASW_G03_2016-07-29",
             "ASW_E04_2016-07-29", "ASW_E11_2016-07-29",
             "ASW_C02_2016-07-29", "Chlorate-PMCs-1-G05_2018-07-01",
             "Chlorate-PMCs-1-E06_2018-07-01", "MK886-PMCs-2-A05_2018-07-01",
             "MK886-PMCs-2-F10_2018-07-01", "MK886-PMCs-3-E05_2018-07-01",
             "MK886-PMCs-3-D01_2018-07-01", "MK886-PMCs-3-B11_2018-07-01",
             "MK886-PMCs-3-G10_2018-07-01", "MK886-PMCs-3-G05_2018-07-01"]

subset_genes = ['SPU_022598', 'SPU_004767', 'SPU_000826', 'SPU_003102',
                'SPU_028395', 'SPU_023386']

def create_annotated_df(expr_file, gene_file, cell_file, filter_cells=None,
                        transpose_expr=True):
    """Create an AnnData object to hold expression, cell, and gene data.
    
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

def filter_cells(anno_data, cell_col, values):
    """
    Filter an annotated data frame to  cells only.
    
    Parameters
    ----------
    anno_data : sc.AnnData
        Annotated data frame with cell annotations denoting treatment.
    cell_col : string
        Name of column in `sc.AnnData.obs` containing features to filter on.
    values : list, string
        list of strings denoting values of interest.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with samples filtered cells.
    """
    keep_cells = [x for x in anno_data.obs.index if\
                  anno_data.obs.loc[x, cell_col] in values]
    return anno_data[keep_cells, ]


def filter_genes(anno_data, gene_col, values):
    """
    Filter an annotated data frame to  cells only.
    
    Parameters
    ----------
    anno_data : sc.AnnData
        Annotated data frame with cell annotations denoting treatment.
    gene_col : string
        Name of column in `sc.AnnData.var` containing features to filter on.
    values : list, string
        list of strings denoting values of interest.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with samples filtered genes.
    """
    keep_genes = [x for x in anno_data.var.index if\
                  anno_data.var.loc[x, gene_col] in values]
    return anno_data[:, keep_genes]


def dispersion(array):
    """[summary]
    
    Parameters
    ----------
    array : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """

    quartiles = np.percentile(array, [25, 75])
    return (quartiles[1] - quartiles[0]) / (quartiles[1] + quartiles[0])

def variable_genes(anno_df, percentile=0.75, ignore_zeros=True):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    percentile : float, optional
        [description] (the default is 0.75, which [default_description])
    ignore_zeros : bool, optional
        [description] (the default is True, which [default_description])
    
    Raises
    ------
    ValueError
        [description]
    
    Returns
    -------
    [type]
        [description]
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


def compare_clusters(anno_df, comparisons=None, cluster_col='louvain',
                     compare_col='treatment', merge_cols=None):
    """
    Compare cell make-up between identified clusters.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe containing cluster membership for each observation
            and a second categorical variable to compare frequency against.
    comparisons : list, optional
        Post-hoc pairwise comparisons to perform. Elements of the list should be
            length 2 containers containing categories to compare (e.g.
            ('Control', 'Expermintanl')). Frequencies between conditions
            will be compared in each cluster using a pairwise G-Test (the
            default is None, which results in no post-hoc analysis). 
    cluster_col : str, optional
        Name of the column in `anno_df.obs` containing observation cluster/group
            assignment (the default is 'louvain'.)
    compare_col : str, optional
        Name of the column in `anno_df.obs` partitioning observations into
            groups of interest (the default is 'treatment', which compares
            distribution of treatments within clusters.)
    merge_cols : dict, optional
        Name of columns to merge together. Keys should point to list of column
            names to merge. The merged column name will be set the key (e.g.
            {'Control': ['ASW', 'DMSO']} will merge 'ASW' and 'DMSO' counts to 
            a single 'Control' column). Default is None, and no merging is
            performed.
    
    Returns
    -------
    [type]
        [description]
    """
    # check for cluster column
    if cluster_col not in anno_df.obs.columns:
        raise ValueError('No {} column in observation data.'.format(
                         cluster_col))

    # check for comparison column
    if compare_col not in anno_df.obs.columns:
        raise ValueError('No {} column in observations data.'.format(
                         compare_col))
    
    #TODO check for comparisons in comparison column

    count_table = pd.crosstab(anno_df.obs[cluster_col],
                              anno_df.obs[compare_col])
    if merge_cols is not None and isinstance(merge_cols, dict):
        for key, columns in merge_cols.items():
            try:
                count_table[key] = count_table[columns].sum(axis=1)
            except KeyError:
                raise('Unknown columns: []'.format(columns))
            count_table.drop(columns, axis=1, inplace=True)
        count_table = count_table[sorted(count_table.columns.values)]
        
    
    # calculate probabilities for comparison values
    probabilities = np.sum(count_table, axis=0) / np.sum(count_table.values)

    # calculate total number of cells per cluster in comparison
    group_counts = np.sum(count_table, axis=1)
    
    # matrix multipy cluster counts and condition probabilities to get
    # expected counts
    expectation = group_counts.values.reshape((count_table.shape[0], 1))\
                  @ probabilities.values.reshape((1, count_table.shape[1]))
    expectation = pd.DataFrame(data=expectation, index=count_table.index,
                               columns=count_table.columns)

    # perform tests of independence between treatments and 
    results = stats.power_divergence(count_table, expectation,
                                     axis=1, lambda_='log-likelihood')

    out = {'gtest': {'observed': count_table, 'expected': expectation,
                     'pvals': results.pvalue,
                     'pvals.adj': results.pvalue * count_table.shape[0]}}

    if comparisons is not None: # perform odds-ratio/fisher exact tests
        fisher_out = {}
        for each in comparisons:
            if len(each) != 2:
                msg = ('Comparisons must be pairwise. Received'
                       ' {} groups: {}'.format(len(each), each))
                raise ValueError(msg)
            
            pairwise = count_table[list(each)]
            out_key = '-'.join(each)
            fisher_out[out_key] = {'odds': np.ones(count_table.shape[0]),
                                   'pvals': np.ones(count_table.shape[0]),
                                   'pvals.adj': np.ones(count_table.shape[0]),
                                   'cluster': count_table.index.values}
            for i, cluster in enumerate(pairwise.index.values):

                # create a 2 x 2 contigency table between treatment comparisons
                # and cluster membership. Counts are pairwise between treatments
                # and cluster X membership vs. not X
                test_cluster = pairwise.loc[cluster, :]
                other_clusters = [x for x in pairwise.index if x != cluster]
                not_cluster = pairwise.loc[other_clusters, :].sum()
                contingency = pd.concat((test_cluster, not_cluster), axis=1)

                # perform fisher exact's test
                odds, pval = stats.fisher_exact(contingency.values)
                fisher_out[out_key]['odds'][i] = odds
                fisher_out[out_key]['pvals'][i] = pval
                fisher_out[out_key]['pvals.adj'][i] = pval * len(comparisons)\
                                                      * count_table.shape[0]
            out['fisher'] = fisher_out
            
    return out

def plot_log_odds(fisher_out):
    
    comparisons = list(fisher_out.keys())
    sorted_comparisons = sorted(comparisons)
    plot_data = pd.DataFrame(fisher_out[comparisons[0]])
    plot_data['Comparison'] = comparisons[0]
    for key in comparisons[1:]:
        new_df = pd.DataFrame(fisher_out[key])
        new_df['Comparison'] = key
        plot_data = pd.concat((plot_data, new_df), axis=0)

    cyndi = True
    if cyndi:
        cyndi_format(len(comparisons))
        plt.grid(False)
        plt.rc('axes.spines', top=False, bottom=False, right=False,
            left=False)
    figure = sns.barplot(data=plot_data, x='cluster', y='odds',
                         hue='Comparison', hue_order=sorted_comparisons)
    plt.axhline(y=1, xmax=plot_data.shape[0], linestyle='--')

    xdiv = 1 / (2*(len(comparisons)) + 1)
    midpoint = len(comparisons) / 2
    for i in range(plot_data.shape[0]):
        row = plot_data.iloc[i, :]
        if row['pvals.adj'] < 0.05:
            marker='*'
            if row['pvals.adj'] < 0.01:
                marker='**'
            pt = np.where(row['Comparison']==np.array(sorted_comparisons))[0][0]
            text_div = 0
            if pt + 1 > midpoint:
                text_div = (pt) * xdiv
            elif pt + 1 < midpoint:
                text_div = -1 * (pt + 1) * xdiv
            elif pt + 1 == midpoint and len(comparisons) % 2 == 0:
                text_div = -1 * (pt + 1) * xdiv
            print(text_div)
            figure.annotate(marker,
                            xy=(row.name + text_div, row['odds']),
                            horizontalalignment='center')

    plt.xlabel('Odds Ratio')
    plt.ylabel('Cluster')

    








def cluster_cells(anno_df, nn=15, metric='cosine'):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    nn : int, optional
        [description] (the default is 15, which [default_description])
    metric : str, optional
        [description] (the default is 'cosine', which [default_description])
    
    Returns
    -------
    [type]
        [description]
    """

    sc.pp.neighbors(anno_df, n_neighbors=nn, metric=metric, use_rep='X')
    sc.tl.umap(anno_df, min_dist=0.0)
    sc.tl.louvain(anno_df)
    anno_df.obs['umap1'] = anno_df.obsm['X_umap'][:, 0]
    anno_df.obs['umap2'] = anno_df.obsm['X_umap'][:, 1]
    # anno_df.obs.to_csv('../../output/plots/clusters.csv')
    sc.pl.umap(anno_df, color='louvain')
    return anno_df


def rename_clusters(anno_df, cluster_col='louvain'):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    cluster_col : str, optional
        [description] (the default is 'louvain', which [default_description])
    
    Returns
    -------
    [type]
        [description]
    """

    clusters = sorted(anno_df.obs[cluster_col].values.unique())
    letters = string.ascii_uppercase[0:len(clusters)]
    int_to_str = {i:a for i, a in zip(clusters, letters)}
    new_clusters = [int_to_str[i] for i in anno_df.obs[cluster_col]]
    anno_df.obs.drop(cluster_col, axis=1, inplace=True)
    anno_df.obs[cluster_col] = new_clusters
    return anno_df


def project_onto_controls(anno_df, neighbors=15):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    neighbors : int, optional
        [description] (the default is 15, which [default_description])
    
    Returns
    -------
    [type]
        [description]
    """

    control_anno = filter_cells(anno_df, 'treatment', ['ASW', 'DMSO'])

    cell_umap = umap.UMAP(n_neighbors=neighbors, metric='cosine', min_dist=0.0)
    cell_umap = cell_umap.fit(control_anno.X)

    all_cords = cell_umap.transform(anno_df.X)
    k_neighbors = int(anno_df.shape[0] / control_anno.shape[0] * neighbors)
    neighbors = sklearn.neighbors.NearestNeighbors(k_neighbors)
    neighbors.fit(all_cords)
    neighbor_graph = neighbors.kneighbors_graph(all_cords)
    anno_df.obsm['X_umap'] = all_cords
    anno_df.obs['umap1'] = all_cords[:, 0]
    anno_df.obs['umap2'] = all_cords[:, 1]
    sc.tl.louvain(anno_df, adjacency=neighbor_graph)
    anno_df = rename_clusters(anno_df, 'louvain')
    return anno_df




    # plt.scatter(cell_umap.embedding_[:, 0], cell_umap.embedding_[:,1], c='blue')
    # plt.scatter(treatment_cords[:, 0], treatment_cords[:, 1], c='red')


    # umap_df = pd.DataFrame(umap_cords, columns=['umap1', 'umap2'])
    # color_dict = {'ASW': 'blue', 'DMSO': 'red'}
    # umap_df.plot.scatter('umap1', 'umap2', c=[color_dict[x] for x in control_anno.obs['treatment']])
    # plt.show()
