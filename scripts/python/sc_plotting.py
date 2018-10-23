import numpy as np
import pandas as pd
import scanpy.api as sc
import seaborn as sns
from matplotlib import pyplot as plt

import sc_utils


def cyndi_format(n_colors):
    """
    Set matplotlib parameters to cyndi's likings.

    Parameters
    ----------
        n_colors : int
            Number of colors to cycle through.
    Returns
    -------
    None
        Sets global plotting parameters to Cyndi's liking.
    """
    from cycler import cycler
    plt.style.use(['dark_background'])
    asw_hex = '#0a66fc'
    chlorate_hex = '#810f7c'
    dmso_hex = '#8cb8ff'
    mk_hex = '#f24343'
    if n_colors == 2:
        plt.rc('axes', prop_cycle=cycler('color', [chlorate_hex, mk_hex]))
    elif n_colors in [3, 4]:
        treatment_colors = [asw_hex, chlorate_hex, dmso_hex, mk_hex]
        if n_colors == 3:
            treatment_colors = [asw_hex, chlorate_hex, mk_hex]
        plt.rc('axes', prop_cycle=cycler('color', treatment_colors))
    else:
        plt.rc('axes', prop_cycle=cycler('color', ['#008fd5', '#fc4f30',
               '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']))

    plt.rcParams["axes.facecolor"] = (0, 0, 0, 0)
    plt.rcParams["figure.facecolor"] = 'black'
    plt.rcParams["axes.labelcolor"] ='white'
    plt.rcParams["text.color"] = "white"
    plt.rcParams["xtick.color"] = "white"
    plt.grid(which='both', color='#c1c1c1', alpha=0.75)
    plt.rcParams['figure.edgecolor'] = '#000000'
    plt.rcParams['axes.facecolor'] = '#000000'


def plot_expectation_difference(observed, expected, normalize=False,
                                compare_col='treatment', p_values=None):
    """
    Plot deviations between observed counts and expected counts.
    
    Parameters
    ----------
    observed : numpy.ndarray
        A (k x p) contingency table of observed counts, where k is the number
           of outcomes/phenotypes, and p is the number treatments/predictors.
    expected : numpy.ndarray
        A (k x p) contingency table of expected counts, where k is the number
        of outcomes/phenotypes, and p is the number treatments/predictors.
    normalize : bool, optional
        Whether to normalize deviations to percent deviations (e.g.
        (expected - count) / expected). The default is False, which plots count
        differences. 
    compare_col : str, optional
        Column name of treatments/predictors. Both `observed` and `expected`
        should have their columns indexed by this name. Default is 'treatment'.
    p_values : list-like, optional
        List of p-values from a test of independence between row and column
        variables (e.g. chi-square, g-test, etc). Values below 0.05 will be
        noted with '*', while values below 0.01 will be noted with '**'. The
        default is None, which produces no annotations.
    
    Returns
    -------
    matplotlib.Axes
        Plot of deviations from expectations.
    """

    plt.style.use('fivethirtyeight')
    plt.rc('axes.spines', top=False, bottom=False, right=False,
            left=False)
    column_order = sorted(observed.columns.values)
    if 'Control' in column_order:
        column_order = ['Control'] + [x for x in column_order if\
                                      x != 'Control']

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
    figure.legend_.set_frame_on(False)

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


def plot_umap(anno_df, color_col=None, shape_col=None):
    """
    Plot cells along their UMAP dimensions.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe with UMAP projections for each cell contained in
        observation data (sc.obs). Projections should be found in columns
        named `umap1` and `umap2`. 
    color_col : str, optional
        Column in observation data (anno_df.obs) determining marker color.
        Default is None, and all observations will be colored the same. 
    shape_col : str, optional
        Column in observation data (anno_df.obs) to determine marker type.
        Default is None, and all observations will be plotted using the same
        marker.
    
    Returns
    -------
    matplotlib.Axes
        Scatter plot of observations along UMAP axes.
    """
    plt.style.use('fivethirtyeight')
    figure = sns.scatterplot(data=anno_df.obs.sort_values(color_col),
                             x='umap1', y='umap2', hue=color_col,
                             style=shape_col, s=100)

    if color_col is not None and shape_col is not None:
        patches, labels = figure.get_legend_handles_labels()
        figure.legend_.remove()
        color_idx = [i for i, x in enumerate(labels) if x == color_col][0]
        shape_idx = [i for i, x in enumerate(labels) if x == shape_col][0]
        labels[color_idx] = labels[color_idx][0].upper() + labels[color_idx][1:]
        labels[shape_idx] = labels[shape_idx][0].upper() + labels[shape_idx][1:]

        color_order = np.argsort(labels[color_idx + 1:shape_idx])
        color_order = [color_idx] + list(color_order + 1)
        color_legend = plt.legend(
                         [patches[color_idx:shape_idx][x] for x in color_order],
                         [labels[color_idx:shape_idx][x] for x in color_order],
                         loc=0, frameon=False)
        shape_order = np.argsort(labels[shape_idx + 1:])
        shape_order = [shape_idx] + list(shape_order + 1)
        figure.legend(
                     [patches[shape_idx:][shape_order][x] for x in color_order],
                     [labels[shape_idx:][shape_order][x] for x in color_order],
                     loc=4, frameon=False)
        plt.gca().add_artist(color_legend)

    elif color_col is not None or shape_col is not None:
        patches, labels = figure.get_legend_handles_labels()
        labels[0] = labels[0][0].upper() + labels[0][1:]
        figure.legend_.remove()
        label_order = np.argsort(labels[1:])
        label_order = [0] + list(label_order + 1)
        figure.legend([patches[x] for x in label_order],
                      [labels[x] for x in label_order], loc=0,
                      frameon=False)
    plt.show()
    plt.cla()


def label_kde(x, color, label):
    """
    Label a ridge plot axis.
    
    Parameters
    ----------
    x : numpy.ndarray
        Values to be plotted.
    color : str
        Color to apply to distribution.
    label : str
        Label to apply to image. 

    Returns
    -------
    matplotlib.Axes
        Axes with ridge plot.
    """
    ax = plt.gca()
    ax.text(0, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)


def ridge_plot(anno_df, gene_name, cluster_col='louvain', mask_zeros=True):
    """
    Create ridge plot of gene expression values striated by group/cluster. 
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe of single-cell expression values.
    gene_name : str
        Gene index id found in `anno_df.var`.
    cluster_col : str, optional
        Column in `anno_df.obs` to group observations by. The default is
        'louvain', which groups cells by louvain clustering results.
    mask_zeros : bool, optional
        Whether to consider zeros during kernel-density estimation. The default
        is True, which removes zeros before kernal estimation.
    
    Returns
    -------
    matplotlib.Axes
        Plot of gene expression striated by group membership.
    """

    # isolate top gene for each plotting
    gene_df = anno_df[:, gene_name]
    gene_df = pd.DataFrame(data=gene_df.X,
                            index=gene_df.obs.index,
                            columns=gene_df.var.index)

    gene_df.columns = ['x']
    if mask_zeros:
        print("Creating ridge plots with masked zero counts. "
              "Set `mask_zeros=False` to keep zeros.")
        gene_df['x'] = gene_df['x'].where(gene_df['x']> 0)
    gene_df[cluster_col] = pd.Series(anno_df.obs[cluster_col], dtype='category')

    # five thirty eight colors
    colors = ['#008fd5', '#fc4f30', '#e5ae38', '#6d904f', '#8b8b8b', '#810f7c']

    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    facet = sns.FacetGrid(data=gene_df, row=cluster_col, hue=cluster_col,
                          palette=colors, aspect=5, height=1,
                          hue_order=sorted(anno_df.obs[cluster_col].unique()))
    facet.map(sns.kdeplot, 'x', clip_on=False, shade=True, alpha=1, lw=1.5,
              bw=.2)
    facet.map(sns.kdeplot, 'x', clip_on=False, color="white", lw=2, bw=.2)
    facet.map(plt.axhline, y=0, lw=2, clip_on=False)

    facet.map(label_kde, 'x')

    # Set the subplots to overlap
    facet.fig.subplots_adjust(hspace=-0.25)

    # Remove axes details that don't play well with overlap
    facet.set_titles("")
    facet.set(yticks=[])
    facet.despine(bottom=True, left=True)
    plt.xlabel('{} expression ($log_2(read counts))$'.format(
               sc_utils.get_gene_identifier(gene_name, anno_df.var)))
    return facet


def plot_de_genes(anno_df, de_results, qthresh=0.001, zero_thresh=0.25,
                  comparison='louvain', n_genes=10):
    """
    Plot top differentially expressed genes.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe with single-cell expression data. 
    de_results : pd.DataFrame
        Dataframe containing differential expression results. Assumed to have a
        'qval' column containing qvalues.
    qthresh : float, optional
        Threshold for 'qval' significance. Default is 0.001.
    zero_thresh : float, optional
        Threshold for zero drop out. Only consider genes with drop out
        percentages above the value. Default is 0.025.
    comparison : str, optional
        Column containing group ids used during differentialy expression. The
        default is 'louvain'.
    n_genes : int, optional
        Number of genes to plot. Default is 10.
    
    Returns
    -------
    None
    """
    perc_zero = np.array([sum(x == 0) / len(x) for x in anno_df.X.T])
    high_genes = anno_df.var.index.values[np.where(perc_zero > zero_thresh)[0]]
    sig_genes = de_results.index[de_results['qval'] < qthresh]
    keep_genes = list(set(high_genes).intersection(sig_genes))
    filtered_de = de_results.filter(keep_genes, axis=0)
    filtered_de.sort_values('qval', inplace=True)
    filtered_anno = anno_df[:, keep_genes]
    filtered_de['gene.name'] = [get_gene_identifier(x, filtered_anno.var)\
                                for x in filtered_de.index.values]
    for gene in filtered_de.index.values[0:n_genes]:
        ridge_plot(filtered_anno, gene, comparison)
        plt.show()
        plt.cla()


def plot_ranked_genes(anno_df, cluster_col='louvain', n_genes=10):
    """
    Plot top-ranked genes in 1-to-all comparisons.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe with ranked genes.
    cluster_col : str, optional
        Column grouping cells. Default is 'louvain'.
    n_genes : int, optional
        Number of genes to grab from each cluster. Default is 10.
    
    Returns
    -------
    None
    """

    ranked_genes = anno_df.uns['rank_genes_groups']
    clusters = sorted(anno_df.obs[cluster_col].values.unique())
    for each in clusters:
        for i in range(n_genes):
            ridge_plot(anno_df, ranked_genes['names'][each][i])


def plot_log_odds(fisher_out):
    """
    Plot log odds between
    
    Parameters
    ----------
    fisher_out : dictionary
        Dictionary containing output from a pairwise Fisher's exact test.
        Generally output from `sc_clustering.compare_clusters()`.
        
        key-value pairs:
            'odds': numpy.array of odds ratios between comparisons ordered by
                    cluster.
            'pvals': numpy.array of calculated p-values orderd by cluster.
            'pvals.adj': numpy.array of adjusted p-values by bonferonni ordered
                         by cluster.
            'cluster': id denoting which cluster comparisons were made in. 

    Returns
    -------
    matplotlib.Axes
        Boxplot of odd-ratios with significant associations denoted with '*'
        and '**' for p-values < 0.05 and p < 0.01, respectively.
    """
    comparisons = list(fisher_out.keys())
    sorted_comparisons = sorted(comparisons)
    plot_data = pd.DataFrame(fisher_out[comparisons[0]])
    plot_data['Comparison'] = comparisons[0]
    for key in comparisons[1:]:
        new_df = pd.DataFrame(fisher_out[key])
        new_df['Comparison'] = key
        plot_data = pd.concat((plot_data, new_df), axis=0)
    figure = sns.barplot(data=plot_data, x='cluster', y='odds',
                         hue='Comparison', hue_order=sorted_comparisons)
    plt.axhline(y=1, xmax=plot_data.shape[0], linestyle='--')

    xdiv = 1 / (2 * (len(comparisons)) + 1)
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
            figure.annotate(marker, xy=(row.name + text_div, row['odds']),
                            horizontalalignment='center')

    plt.xlabel('Odds Ratio')
    plt.ylabel('Cluster')
    return figure


def plot_mixture(X, gmm, gene_name=None):
    """
    Plot the noise and signal mixture model.
    
    Parameters
    ----------
    X : numpy.array
        Array of single-cell gene expression profile for a single gene.
    gmm : pomegranate.GeneralMixtureModel
        General mixture model modelling noise and signal portions of
        single-cell expression profile.
    gene_name : str, optional
        Name of gene being modeled. Name will appear in plot title. Default is
        None, and no gene name will appear.
    
    Returns
    -------
    matplotlib.Axes
        Plot of noise-signal mixture model.
    """
    if len(X.shape) == 1:
        X = X.reshape(-1, 1)
    sns.set_palette('deep')
    colors = sns.color_palette()
    space = np.arange(0, np.max(X), 0.1)
    sns.distplot(a=X, bins=30, kde=False, rug=False, norm_hist=True,
                 color=colors[4])
    noise_lambda = gmm.distributions[0].parameters[0]
    noise_model = 'Poisson($\lambda={:0.2f}$)'.format(noise_lambda)
    mu, sig = gmm.distributions[1].parameters
    signal_model = "N($\mu={:0.2f}$, $\sigma={:0.2f}$)".format(mu, sig)
    plt.plot(space, gmm.distributions[0].probability(space), linestyle=':',
             label='Noise Model ~ {}'.format(noise_model), color=colors[3],
             linewidth=3)
    plt.plot(space, gmm.distributions[1].probability(space), linestyle=':',
             label='Signal Model ~ {}'.format(signal_model), color=colors[0],
             linewidth=3)
    filtered, threshold = sc_utils.threshold_expression(X, method='mixture')
    plt.axvline(x=threshold, linestyle='--',
                label='Noise Threshold = {}'.format(threshold), color='black')
    plt.plot(space, gmm.probability(space), label='Mixture Model',
             color=colors[4], linewidth=3)
    plt.legend(fontsize=14, loc=0)
    plt.xlabel('$\log_2(read counts)$', fontsize=16)
    plt.ylabel('$p(x)$', fontsize=16)
    title = "Signal-Noise Mixture Model"
    if gene_name != None:
        title += ' for {}'.format(gene_name)
    plt.title(title, loc='left', fontsize=18)


def expression_heatmap(anno_df, genes, cluster_cells=True, cluster_genes=True):
    """[summary]
    
    Parameters
    ----------
    anno_df : [type]
        [description]
    genes : [type]
        [description]
    cluster_cells : bool, optional
        [description] (the default is True, which [default_description])
    cluster_genes : bool, optional
        [description] (the default is True, which [default_description])
    
    """

    plot_data = anno_df[:, genes]
    sns.clustermap(plot_data.X, cmap='inferno', row_cluster=cluster_cells,
                   col_cluster=cluster_genes, standard_scale=1)
    plt.show()


def distribution_plots(X):
    """[summary]
    
    Parameters
    ----------
    X : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
    sns.distplot(X, bins=30, ax=ax1)
    sns.distplot(X, bins=30, ax=ax2, hist_kws={'cumulative': True},
                kde_kws={'cumulative': True})

    return fig


