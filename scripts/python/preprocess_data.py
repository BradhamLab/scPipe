"""
Filter count matrices from scRNAseq data to remove uninformative genes and
cells.

author: Dakota Hawkins
date: August 20, 2018. 
"""

import os
import sys 

import pandas as pd

def filter_count_matrix(count_matrix, n_reads=50000, n_cells=2, cpm_thresh=0.5):
    """
    Remove poorly sequenced cells and genes with low occurrence.

    Filter cells if the total aligned read counts falls below a provided
    threshodl. Filter genes with both low counts-per-million and low occurrence
    among cells.

    Args:
        count_matrix (pd.DataFrame): raw read-count matrix.
        n_reads (int): minimum number of read counts a cell must have to avoid
            filtering. Default is 50,000.
        n_cells (int, optional): minimum number of cells required to exhibit
            the minimum expression level to avoid filtering. Default is 2. 
        cpm_thresh (float, optional): minimum counts-per-million in lowly
            mapped genes. Default is 0.5
    Returns:
        (pd.DataFrame): filtered dataframe.

    References:
        Cell Filtering:
            Rizzetto, S., Eltahla, A. A., Lin, P., Bull, R., Lloyd, A. R., Ho,
            J. W. K., … Luciani, F. (2017). Impact of sequencing depth and read
            length on single cell RNA sequencing data of T cells.
            Scientific Reports, 7(1), 12781.
            https://doi.org/10.1038/s41598-017-12989-x
            https://www.nature.com/articles/s41598-017-12989-x

        Gene Filtering:
            Chen Y, Lun ATL and Smyth GK. From reads to genes to pathways:
            differential expression analysis of RNA-Seq experiments using
            Rsubread and the edgeR quasi-likelihood pipeline
            [version 2; referees: 5 approved]. F1000Research 2016, 5:1438
            (doi: 10.12688/f1000research.8987.2) 
            https://f1000research.com/articles/5-1438/v2
        
    """
    # drop cells with low coverage
    cell_counts = count_matrix.sum()
    bad_cells = cell_counts.index.values[cell_counts < n_reads]
    count_matrix.drop(bad_cells, axis=1, inplace=True)

    # drop genes with low expression and low occurrence
    cpm = count_matrix.apply(lambda x: x / cell_counts[x.name] * 10**6,
                              axis=0)
    low_genes = cpm.apply(lambda x: sum(x > cpm_thresh) < n_cells, axis=1)
    low_genes = low_genes.index.values[low_genes]
    return count_matrix.drop(labels=low_genes, axis=0)


def filter_metadata(metadata, cells):
    """
    Remove filtered cell from metadata dataframe.
    
    Args:
        metadata (pd.DataFrame): cell metadata.
        cells (list-like): list of cells to keep.
    
    Returns:
        (pd.DataFrame): filtered metadata.
    """
    
    return metadata.loc[cells]

def main(count_matrix, metadata, n_reads=50000, n_cells=2, cpm_thresh=0.5):
    """Filter count and meta data."""
    filtered_count = filter_count_matrix(count_matrix, n_reads, n_cells,
                                         cpm_thresh)
    
    filtered_meta = filter_metadata(metadata, filtered_count.columns.values)

    return filtered_count, filtered_meta


if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False

    if snakemake_exists:
        cmat = pd.read_csv(snakemake.input['cmat'], index_col=0)
        tpm = pd.read_csv(snakemake.input['tpm'], index_col=0)
        meta = pd.read_csv(snakemake.input['meta'], index_col=0)
        f_cmat, f_meta = main(cmat, meta, snakemake.params['reads'],
                              snakemake.params['cells'],
                              snakemake.params['cpm'])
        f_tpm = tpm.loc[f_cmat.index.values, f_cmat.columns.values]
        f_cmat.to_csv(snakemake.output['cmat'])
        f_tpm.to_csv(snakemake.output['tpm'])
        f_meta.to_csv(snakemake.output['meta'])