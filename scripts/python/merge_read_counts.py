"""
Merge sample count files into a count matrix with genes as rows and samples as
columns.

Author: Dakota Hawkins
Date: August 16, 2018
"""

import os

import pandas as pd
import numpy as np

def merge_counts(count_dir):
    """
    Merge featureCounts output from several sample files into a
    single dataframe.

    Merge featureCounts quantification across samples. Returns raw feature
    counts across datasets, as well as transcripts per million (TPM) counts.

    Args:
        count_dir (string): path to directory containing HTSeq output. Files
        are expected to follow a {sample.name}.txt naming scheme.

    Returns:
        pd.DataFrame: Pandas dataframes containing raw read counts
    """
    data_frames = []
    for x in os.listdir(count_dir):
        sample = os.path.basename(os.path.splitext(x)[0])
        sample_df = pd.read_csv(os.path.join(count_dir, x),
                                index_col=0, sep='\t',
                                header=0,
                                names=['gene_id', sample],
                                skiprows=[0])
        data_frames.append(sample_df[sample])

    return pd.concat([each for each in data_frames], axis=1)


def remove_zero_genes(count_df):
    """Remove genes with no reads mapping to them across samples."""
    count_df[count_df == 0] = np.NaN
    filtered = count_df.dropna(how='all')
    return filtered.fillna(value=0)


if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None

    if snakemake is not None:
        counts = merge_counts(snakemake.params['dir'])
        filtered_counts = remove_zero_genes(counts)
        counts.to_csv(snakemake.output['count'])
