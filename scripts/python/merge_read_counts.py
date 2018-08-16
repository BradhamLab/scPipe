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

    Args:
        count_dir (string): path to directory containing HTSeq output. Files
        are expected to follow a {sample.name}.txt naming scheme.

    Returns:
     (pd.DataFrame): pandas dataframe of read counts mapped to each gene.
    """
    data_frames = []
    for x in os.listdir(count_dir):
        data_frames.append(pd.read_csv(os.path.join(count_dir, x),
                                       index_col=0, sep='\t',
                                       names=[os.path.splitext(x)[0]],
                                       skiprows=[0, 1]))
    count_df = pd.concat(data_frames, axis=1)
    return count_df


def remove_zero_genes(count_df):
    """
    Remove genes with no reads mapping to them across samples.
    """
    count_df[count_df == 0] = np.NaN
    filtered = count_df.dropna(how='all')
    return filtered.fillna(value=0)


if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False

    if snakemake_exists:
        count_df = merge_counts(snakemake.params['dir'])
        filtered_df = remove_zero_genes(count_df)
        filtered_df.to_csv(snakemake.output['csv'])
