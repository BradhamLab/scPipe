import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import scanpy.api as sc
import umap
import louvain
from scipy import spatial

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
    