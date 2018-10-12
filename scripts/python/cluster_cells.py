import string

import igraph
import louvain
import matplotlib.pyplot as plt
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

# Define and use a simple function to label the plot in axes coordinates
def label_kde(x, color, label):
    """[summary]
    
    Parameters
    ----------
    x : [type]
        [description]
    color : [type]
        [description]
    label : [type]
        [description]
    
    """

    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


# TODO document
def get_gene_identifier(scaffold, gene_df):
    """[summary]
    
    Parameters
    ----------
    scaffold : [type]
        [description]
    gene_df : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """

    if not pd.isnull(gene_df.loc[scaffold, 'UniProt.Name']):
        return gene_df.loc[scaffold, 'UniProt.Name'].split('_')[0]
    elif not pd.isnull(gene_df.loc[scaffold, 'Uniprot.ID']):
        return gene_df.loc[scaffold, 'UniProt.ID']
    elif not pd.isnull(gene_df.loc[scaffold, 'SPU']):
        return gene_df.loc[scaffold, 'SPU']
    else:
        return "No.Name"


#TODO document
def ridge_plot(anno_df, gene_name, cluster_col='louvain', mask_zeros=True):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    gene_name : [type]
        [description]
    cluster_col : str, optional
        [description] (the default is 'louvain', which [default_description])
    mask_zeros : bool, optional
        [description] (the default is True, which [default_description])
    
    """

    # isolate top gene for each plotting
    gene_df = anno_df[:, gene_name]
    gene_df = pd.DataFrame(data=gene_df.X,
                            index=gene_df.obs.index,
                            columns=gene_df.var.index)
    gene_df.columns = ['x']
    if mask_zeros:
        print("Creating ridge plots with masked zero counts. Set `mask_zeros=False` to keep zeros.")
        gene_df['x'] = gene_df['x'].where(gene_df['x']> 0)
    gene_df[cluster_col] = anno_df.obs[cluster_col]

    # five thirty eight colors
    colors = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    facet = sns.FacetGrid(data=gene_df, row=cluster_col, hue=cluster_col,
                          palette=colors, aspect=5, height=1)
    facet.map(sns.kdeplot, 'x', clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    facet.map(sns.kdeplot, 'x', clip_on=False, color="w", lw=2, bw=.2)
    facet.map(plt.axhline, y=0, lw=2, clip_on=False)

    facet.map(label_kde, 'x')

    # Set the subplots to overlap
    facet.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    facet.set_titles("")
    facet.set(yticks=[])
    facet.despine(bottom=True, left=True)
    plt.xlabel('{} expression ($log_2(read counts))$'.format(
        get_gene_identifier(gene_name, anno_df.var)))
    plt.show()
    plt.cla()


# TODO document
def plot_de_genes(anno_df, de_results, comparison='louvain', n_genes=10):
    perc_zero = [sum(x != 0) / len(x) for x in anno_df.X]



#TODO document
def plot_ranked_genes(anno_df, cluster_col='louvain', n_genes=10):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    cluster_col : str, optional
        [description] (the default is 'louvain', which [default_description])
    n_genes : int, optional
        [description] (the default is 10, which [default_description])
    
    """

    ranked_genes = anno_df.uns['rank_genes_groups']
    clusters = sorted(anno_df.obs[cluster_col].values.unique())
    best_genes = set()
    for each in clusters:
        ridge_plot(anno_df, ranked_genes['names'][each][0])
        best_genes = best_genes.union(ranked_genes['names'][each])

    # heatmap_data = pd.DataFrame(data=anno_df[:, list(best_genes)].X,
    #                             index=anno_df.obs.index,
    #                             columns=list(best_genes))
    # colors = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']
    # clus2col = {x: colors[i % len(clusters)] for i, x in enumerate(clusters)}
    # cluster_colors = [clus2col[x] for x in anno_df.obs[cluster_col]]
    # sns.clustermap(data=heatmap_data, row_colors=cluster_colors)


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


# TODO document
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
    cyndi = False
    plt.rc('axes.spines', top=False, bottom=False, right=False,
            left=False)

    column_order = sorted(observed.columns.values)
    if 'Control' in column_order:
        column_order = ['Control'] + [x for x in column_order if\
                                      x != 'Control']
    if cyndi:
        cyndi_format(observed.shape[1])
        difference = observed.apply(lambda x: x / np.sum(x) * 100, 0)
        difference = difference[column_order]
        ylab = 'Percentage Count'

    else:
    plt.style.use('fivethirtyeight')
    if not normalize:
        # raw count difference
        difference = observed - expected
        ylab = 'Count Difference'
    else:
        # calculate percent error 
        difference = (observed - expected) / (expected) * 100
            ylab = 'Deviation from Expectation'
    
    # treat categorical data as strings
    difference.rename(columns=str, inplace=True)
    difference['Cluster'] = difference.index.values
    plot_data = pd.melt(difference, id_vars='Cluster')
    figure = sns.barplot(x='Cluster', y='value', hue=compare_col,
                         data=plot_data, hue_order=column_order)

    # fix axes
    if normalize:
        ylocs, ylabs = plt.yticks()
        new_labs = ['{:.0f}%'.format(each) for each in ylocs]
        plt.yticks(ylocs, new_labs)
    plt.ylabel(ylab)

    # fix legend
    figure.legend_.set_title(compare_col[0].upper() + compare_col[1:])
    # figure.legend_.set_frame_on(False)

    # add denetors for statistically significant deviations
    
    if p_values is not None:
        ydiv = max((plot_data['value'].max() - plot_data['value'].min()) * 0.01,
                    0)
        for i, p in enumerate(p_values):
            if p < 0.05:
                marker='*'
                if p < 0.01:
                    marker='**'
                figure.annotate(marker,
                                xy=(i, np.max(difference.iloc[i, :-1])),
                                horizontalalignment='center')
            n = observed.iloc[i, :].sum()
            text_div = -1 * np.sign(difference.iloc[i, :-1][-1]) * ydiv
            vert_align = 'top'
            if  text_div > 0:
                vert_align = 'bottom'
            figure.annotate('(n={})'.format(n),
                             xy=(i, 0),
                             xycoords='data',
                             xytext=(i + 0.5, text_div),
                             textcoords='data',
                             horizontalalignment='right',
                             verticalalignment=vert_align,
                             size='x-small')
    
    return figure

def plot_umap(anno_df, color_col, shape_col):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    color_col : [type]
        [description]
    shape_col : [type]
        [description]
    
    """

    plt.style.use('fivethirtyeight')
    figure = sns.scatterplot(data=clustered.obs, x='umap1', y='umap2',
                             hue=color_col, style=shape_col, s=100)
    patches, labels = figure.get_legend_handles_labels()
    figure.legend_.remove()

    color_idx = [i for i, x in enumerate(labels) if x == color_col][0]
    shape_idx = [i for i, x in enumerate(labels) if x == shape_col][0]
    labels[color_idx] = labels[color_idx][0].upper() + labels[color_idx][1:]
    labels[shape_idx] = labels[shape_idx][0].upper() + labels[shape_idx][1:]

    label_order = np.argsort(labels[color_idx:shape_idx])
    color_legend = plt.legend(patches[color_idx:shape_idx],
                              labels[color_idx:shape_idx], loc=0,
                              frameon=False)
    figure.legend(patches[shape_idx:], labels[shape_idx:], loc=4, frameon=False)
    plt.gca().add_artist(color_legend)
    plt.show()
    plt.cla()


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
