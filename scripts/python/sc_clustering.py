import string

import numpy as np
import pandas as pd
import scanpy.api as sc
import sklearn
import umap
from scipy import spatial, stats

import sc_utils

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
    dictionary
        Dictionary containing comparison results between clusters including
        G-tests of independence and, if a `comparisons` argument was provided,
        pairwise Fisher exact test results.

        key-value pairs:
            'gtest': dictionary of Gtest results
                key-value pairs:
                    'observed': matrix of observed counts.
                    'expected': matrix of expected counts.
                    'pvals': pvalues for a GTest of independence following
                        index order of 'observed' and 'expected' tables.
                    'pvals.adj': adjusted p-values using the bonferonni
                        correction.
            'fisher': dictionary of pairwise comparison results. Results are
                contained in a dictionary keyed by comparison groups separated
                by a hyphen (e.g. if ['X', 'Y'] was provided as a comparison,
                results of the comparison would be keyed 'X-Y').
                
                key-value pairs:
                    'odds': numpy.array of odds ratios between comparisons
                        ordered by cluster.
                    'pvals': numpy.array of calculated p-values orderd by
                        cluster.
                    'pvals.adj': numpy.array of adjusted p-values by bonferonni
                        ordered by cluster.
                    'cluster': id denoting which cluster comparisons. 
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
                raise('Unknown columns: {}'.format(columns))
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

    
def cluster_cells(anno_df, nn=15, metric='cosine'):
    """
    Cluster cells using Louvain modularity on UMAP projections of cells.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe containing expression data
    nn : int, optional
        Number of neighbors to use in UMAP projections. Default is 15.
    metric : str, optional
        Distance metric to use. Default is 'cosine'.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with clustering results.
    """

    sc.pp.neighbors(anno_df, n_neighbors=nn, metric=metric)
    sc.tl.umap(anno_df, min_dist=0.0)
    sc.tl.louvain(anno_df)
    anno_df.obs['umap1'] = anno_df.obsm['X_umap'][:, 0]
    anno_df.obs['umap2'] = anno_df.obsm['X_umap'][:, 1]
    # anno_df.obs.to_csv('../../output/plots/clusters.csv')
    sc.pl.umap(anno_df, color='louvain')
    return anno_df


def rename_clusters(anno_df, cluster_col='louvain'):
    """
    Rename cluster ids to single letter ids.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated data frame with clustered cells.
    cluster_col : str, optional
        Column in observations data (`anno_df.obs`) containing cluster
        assignments.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with re-named cluster ids.
    """

    clusters = sorted(anno_df.obs[cluster_col].values.unique())
    letters = string.ascii_uppercase[0:len(clusters)]
    int_to_str = {i:a for i, a in zip(clusters, letters)}
    new_clusters = [int_to_str[i] for i in anno_df.obs[cluster_col]]
    anno_df.obs.drop(cluster_col, axis=1, inplace=True)
    anno_df.obs[cluster_col] = new_clusters
    return anno_df


def project_onto_controls(anno_df, neighbors=15):
    """
    Cluster cells by first learning a projection using control cells. 

    Cluster cells using Louvain modularity on UMAP projected cells. However, 
    UMAP projection is first learned using only control cells. All cells
    are then transformed into the learned space, and clustered using Louvain
    modularity. 
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of expression data.
    neighbors : int, optional
        Number of neighbors to use when learning UMAP transformation. Default is
        15.
    
    Returns
    -------
    sc.AnnData
        Annotated dataframe with clustered cells. Clusters found in 'louvain'
        column in observation dataframe.
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