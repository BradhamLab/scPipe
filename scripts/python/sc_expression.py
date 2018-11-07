import numpy as np
from scipy import stats
from anndata.base import ArrayView
from skimage.filters import threshold_otsu
import re
import pandas as pd

import sc_classes
import sc_utils
import sc_plotting

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
                               compare_zeros=True, plot=True):
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
    remove_noise = True
    if remove_noise:
        X = threshold_expression(anno_df[:, gene].X, method='mixture')[0]
        anno_df[:, gene].X = X

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
        g1_zeros = group_1.shape[0] - g1_not
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
    out_dict['group1.avg'] = float(np.mean(non_zero_g1))
    out_dict['group2.avg'] = float(np.mean(non_zero_g2))
    out_dict['group1.std'] = float(np.std(non_zero_g1))
    out_dict['group2.std'] = float(np.std(non_zero_g2))
    out_dict['gene'] = gene
    out_dict.update(anno_df.var.loc[gene, :].to_dict())
    if plot:
        import seaborn as sns
        plot_df = sc_utils.filter_cells(anno_df, group_col,
                lambda x: any([True for y in regex if y.search(x) is not None]))
        exprs = np.array(plot_df[:, gene].X)
        expr_col = sc_utils.get_gene_identifier(gene, plot_df.var)
        plot_df.obs[expr_col] = exprs
        group_column = []
        for x in plot_df.obs[group_col]:
            if regex[0].search(x) is not None:
                group_column.append(groups[0])
            else:
                group_column.append(groups[1])
        plot_df.obs[group_col] = group_column
        # format count data for plotting
        count_df = pd.DataFrame(data=counts, columns=groups,
                                index=['Zero', 'Not.Zero'])
        count_df.index.name = 'Zero.Status'
        count_df.columns.name= 'Treatment'

        # calculate probabilities for comparison values
        probabilities = np.sum(count_df, axis=0) / np.sum(count_df.values)

        # calculate total number of cells per cluster in comparison
        group_counts = np.sum(count_df, axis=1)
        
        # matrix multipy cluster counts and condition probabilities to get
        # expected counts
        expectation = group_counts.values.reshape((count_df.shape[0], 1))\
                    @ probabilities.values.reshape((1, count_df.shape[1]))
        expectation = pd.DataFrame(data=expectation, index=count_df.index,
                                   columns=count_df.columns)
        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        # plot expectation deviation
        fisher_plot = sc_plotting.plot_expectation_difference(count_df,
                                      expectation, normalize=False,
                                      ax=axes[0][1])
        # plot expression along umap coordinates
        umap_exprs_plot = sc_plotting.plot_umap(plot_df, shape_col=group_col,
                                                gene=gene, ax=axes[0][0])

        # plot expression between treatment
        zero_expr_plot = sns.boxplot(data=plot_df.obs, x=group_col, y=expr_col,
                                     showfliers=False, ax=axes[1][0])
        zero_expr_plot = sns.swarmplot(data=plot_df.obs, x=group_col,
                                       y=expr_col, ax=zero_expr_plot,
                                       color="0.25")
        non_zero = plot_df.obs[plot_df.obs[expr_col] != 0]
        non_zero_plot = sns.boxplot(data=non_zero, x=group_col, y=expr_col,
                                    showfliers=False, ax=axes[1][1])
        non_zero_plot = sns.swarmplot(data=non_zero, x=group_col, y=expr_col, 
                                      ax=non_zero_plot, color="0.25")
        # file_name = "../../output/plots/ExprPlots/" + expr_col + '.png'
        fig.savefig(file_name)
        plt.close()
        
    return out_dict


if __name__ == "__main__":
    import sc_clustering
    anno_df = sc_utils.create_annotated_df(expr_file=sc_clustering.expr_file,
                                           gene_file=sc_clustering.gene_file,
                                           cell_file=sc_clustering.cell_file,
                                           filter_cells=sc_clustering.bad_cells)
    dan_genes = pd.read_csv('../../files/dan_genes.csv')
    dan_spus = dan_genes['Dan.SPU.updates'].values
    dan_df = sc_utils.filter_genes(anno_df, 'SPU', lambda x: x in dan_spus)
    genes = dan_df.var.index.values