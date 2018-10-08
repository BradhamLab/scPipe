import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy.api as sc
import umap
import louvain
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

def create_annotated_df(expr_file, gene_file, cell_file, filter_cells=None):
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
    
    Returns
    -------
    sc.AnnData
        Annotated data frame containing expression data, cell data, and gene
        data.
    """

    # read in data
    expr_data = pd.read_csv(expr_file, index_col=0)
    gene_data = pd.read_csv(gene_file, index_col=0)
    gene_data['scaffold'] = [x.replace('model', 'TU') for x in gene_data.index]
    gene_data = gene_data.set_index('scaffold')  # this is filling with NaN

    cell_data = pd.read_csv(cell_file, index_col=0)

    # remove non-PMC cells
    if filter_cells is not None:
        keep_cells = [x for x in cell_data.index if x not in filter_cells]
        cell_data = cell_data.filter(keep_cells, axis=0)
        expr_data = expr_data.filter(keep_cells, axis=1)

    # remove genes that are not in both annotation and and expression datasets 
    keep_genes = list(set(expr_data.index).intersection(gene_data.index))
    gene_data = gene_data.filter(keep_genes, axis=0)
    expr_data = expr_data.filter(keep_genes, axis=0)

    # create anno df formatted cells x genes 
    anno_df = sc.AnnData(expr_data.T.values, obs=cell_data, var=gene_data)
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


def compare_clusters(anno_df, comparisons=None, cluster_col='louvain',
                     compare_col='treatment'):
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
    
    # calculate probabilities for comparison values
    probabilites = np.sum(count_table, axis=0) / np.sum(count_table.values)

    # calculate total number of cells per cluster in comparison
    group_counts = np.sum(count_table, axis=1)
    
    # matrix multipy cluster counts and condition probabilities to get
    # expected counts
    expectation = group_counts.values.reshape((count_table.shape[0], 1))\
                  @ probabilites.values.reshape((1, count_table.shape[1]))
    expectation = pd.DataFrame(data=expectation, index=count_table.index,
                               columns=count_table.columns)

    # TODO do comparison off of original p-values
    # TODO check pairwise implementation 
    if comparisons is not None:
        out = {}
        for each in comparisons:
            # perform g-test. TODO: look into degrees of freedom
            results = stats.power_divergence(count_table[list(each)],
                                             expectation[list(each)],
                                             axis=1, lambda_='log-likelihood')

            adj_pvals = results.pvalue * len(comparisons)\
                                       * count_table.shape[0]

            out['-'.join(list(each))] = {'observed': count_table[list(each)],
                                         'expected': expectation[list(each)],
                                         'pvals': results.pvalue,
                                         'pvals.adj': adj_pvals}
    else:
        results = stats.power_divergence(count_table, expectation,
                                            axis=1, lambda_='log-likelihood')
        adj_pvals = results.pvalue * count_table.shape[0]

        out = {'observed': count_table, 'expected': expectation,
                'pvals': results.pvalue, 'pvals.adj': adj_pvals}

    return out


def plot_expectation_difference(observed, expected, normalize=False,
                                compare_col='treatment', p_values=None):
    """
    Plot 
    
    Parameters
    ----------
    observed : [type]
        [description]
    expected : [type]
        [description]
    normalize : bool, optional
        [description] (the default is False, which [default_description])
    compare_col : str, optional
        [description] (the default is 'treatment', which [default_description])
    p_values : [type], optional
        [description] (the default is None, which [default_description])
    
    Returns
    -------
    [type]
        [description]
    """

    plt.style.use('fivethirtyeight')
    if not normalize:
        # raw difference
        difference = observed - expected
        ylab = 'Count Difference'
    else:
        # calculate percent error 
        difference = (observed - expected) / (expected) * 100
        ylab = 'Percent Error'
        
    difference['Cluster'] = difference.index.values
    plot_data = pd.melt(difference, id_vars='Cluster')
    figure = sns.barplot(x='Cluster', y='value', hue=compare_col,
                         data=plot_data)

    # fix axes
    if normalize:
        ylocs, ylabs = plt.yticks()
        new_labs = ['{:.0f}%'.format(each) for each in ylocs]
        plt.yticks(ylocs, new_labs)
    plt.ylabel(ylab)

    # fix legend
    figure.legend_.set_title(compare_col[0].upper() + compare_col[1:])
    figure.legend_.set_frame_on(False)

    # add denetors for statistically significant deviations
    if p_values is not None:
        for i, p in enumerate(p_values):
            if p < 0.5:
                marker='*'
                if p < 0.01:
                    marker='**'
                figure.annotate(marker,
                                xy=(difference.index.values[i],
                                    np.max(difference.iloc[i, :-1])),
                                horizontalalignment='center')
    
    return figure

def cluster_cells(anno_df, nn=15, metric='cosine'):
    sc.pp.neighbors(anno_df, n_neighbors=nn, metric=metric)
    sc.tl.umap(anno_df, min_dist=0.0)
    sc.tl.louvain(anno_df)
    anno_df.obs['umap1'] = anno_df.obsm['X_umap'][:, 0]
    anno_df.obs['umap2'] = anno_df.obsm['X_umap'][:, 1]
    # anno_df.obs.to_csv('../../output/plots/clusters.csv')
    sc.pl.umap(anno_df, color='louvain')
    return anno_df
    # treatments = anno_df.obs['treatment'].unique()
    # clusters = anno_df.obs['louvain'].unique()

    # anno_df.obs.plot.scatter('umap1', 'umap2', anno_df.obs['louvain'],
    #                          marker=anno_df.obs['treatment'])
    # cell_umap = umap.UMAP(metric='cosine')
    # cell_umap = cell_umap.fit(control_anno.X)
    # treatment_cords = cell_umap.transform(treatment_df.X)
    # embeddings = np.vstack((cell_umap.embedding_, treatment_cords))
    # cells = np.hstack((control_anno.obs.index.values,
    #                    treatment_df.obs.index.values))
    # combined_expr = np.vstack((control_anno.X, treatment_df.X))
    # combined_obs = pd.concat((control_anno.obs, treatment_df.obs), axis=0)
    # combined_obs['umap1'] = embeddings[:, 0]
    # combined_obs['umap2'] = embeddings[:, 1]


    # plt.scatter(cell_umap.embedding_[:, 0], cell_umap.embedding_[:,1], c='blue')
    # plt.scatter(treatment_cords[:, 0], treatment_cords[:, 1], c='red')


    # umap_df = pd.DataFrame(umap_cords, columns=['umap1', 'umap2'])
    # color_dict = {'ASW': 'blue', 'DMSO': 'red'}
    # umap_df.plot.scatter('umap1', 'umap2', c=[color_dict[x] for x in control_anno.obs['treatment']])
    # plt.show()
    