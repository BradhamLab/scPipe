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
     (dict, pd.DataFrame): dictionary of pandas dataframes containing raw read
        counts and TPM normalized values. Keyed by 'count' and 'tpm'. 
    """
    data_frames = []
    for x in os.listdir(count_dir):
        sample = os.path.basename(os.path.splitext(x)[0])
        sample_df = pd.read_csv(os.path.join(count_dir, x),
                                index_col=0, sep='\t',
                                header=0,
                                names=['Length', sample],
                                skiprows=[0])
        # calculate count to feature length ratio
        counts_per_base = sample_df.apply(lambda x: x[1] / x[0], axis=1)

        # calculate sum of ratios / rates for normalization
        total_cpb = counts_per_base.sum()
        if total_cpb != 0:
            tpm = counts_per_base.apply(lambda x: x * 1 / total_cpb * 10**6)
            tpm.rename(sample, inplace=True)
        # total rates is zero -> no transcriptcs, tpm = 0 for all 
        else:
            tpm = sample_df[sample]

        data_frames.append({'count': sample_df[sample], 'tpm': tpm})

    count_matrix = pd.concat([each['count'] for each in data_frames], axis=1)
    tpm_matrix = pd.concat([each['tpm'] for each in data_frames], axis=1)
    return {'count':count_matrix, 'tpm': tpm_matrix}


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
        dataframes = merge_counts(snakemake.params['dir'])
        filtered_counts = remove_zero_genes(dataframes['count'])
        filtered_tpm = dataframes['tpm'].loc[filtered_counts.index,
                                             filtered_counts.columns]
        filtered_counts.to_csv(snakemake.output['count'])
        filtered_tpm.to_csv(snakemake.output['tpm'])
