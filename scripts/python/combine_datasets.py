"""
Combine count and metadata matrices from different single-cell RNAseq runs.

Author: Dakota Hawkins
Date: August 16, 2018
"""

# system imports 
import sys
import os

import itertools

import pandas as pd

def load_datasets(output_dirs):
    """
    Load all dataset to be merged.

    Args:
        output_dirs (list, string): list of scPipe output directories.
    Return:
        (list, dict): list of dictionaries keyed by 'counts' and 'meta' to hold
            count matrices and sample metadata for each run.
    """
    datasets = []
    for each in output_dirs:
        cmatrix = pd.read_csv(os.path.join(each, 'matrix', 'count_matrix.csv'),
                              index_col=0)
        meta = pd.read_csv(os.path.join(each, 'metadata', 'metadata.csv'),
                           index_col=0)
        cmatrix, meta = ensure_unique_batch_and_sample_ids(cmatrix, meta)
        datasets.append({'counts': cmatrix, 'meta': meta})

    return datasets


def ensure_unique_batch_and_sample_ids(cmat, meta):
    """
    Ensure unique sample ids and batch ids between datasets.

    Concatenates `run.id` with index values and batch ids to ensure unique
    identifiers between datasets.

    Args:
        cmat (pd.DataFrame): count matrix.
        meta (pd.DataFrame): sample data.
    Returns:
        (tuple, (pd.DataFrame, pd.DataFrame)): count and sample data matrices
            with unique sample ids and batch ids. (cmat, meta).
    """

    # exit if no 
    if 'run.id' not in meta.columns:
        sys.exit('Error: sample data matrix must contain run identifying\
\\ \\ \\ \\ \\ \\information in a `run.id` column')

    run_id = meta['run.id'].values[0]

    # ensure batch column exists 
    if 'batch' not in meta.columns:
        meta['batch'] = run_id

    # ensure unique batch ids across runs
    meta['batch'] = meta.apply(lambda x: '{}-{}'.format(run_id, x['batch']),
                                  axis=1)

    # ensure unique sample ids in meta data across runs
    meta['new_index'] = meta.apply(lambda x: '{}_{}'.format(x.name, run_id),
                                   axis=1)
    meta.set_index('new_index', inplace=True)

    # ensure unique sample names in count data across runs
    new_columns = cmat.apply(lambda x: '{}_{}'.format(x.name, run_id), axis=0)
    cmat.columns = new_columns

    return(cmat, meta)


def rewrite_originals(cmat_csv, meta_csv):
    """
    Rewrite original data to new files.

    Rewrites original count and meta data to new files. Given the original
    file path {/foo/bar/count_matrix.csv}, the original data will be written to
    {/foo/bar/original_count_matrix.csv}.

    Args:
        cmat_csv (string): path to count matrix file.
        meta_csv (string): path to meta data file.
    Returns:
        None
    """
    for each in [cmat_csv, meta_csv]:
        data = pd.read_csv(each, index_col=0)
        path_split = os.path.split(each)
        new_file = os.path.join(path_split[0], 'original_' + path_split[1])
        data.to_csv(new_file)


def combine_datasets(datasets):
    """
    Combine count and meta data between datasets.

    See `combine_count_data()` and `combine_meta_data()` for more information.

    Args:
        datasets (list, dict): list of dictionaries containing count and meta
            data for each alignment run.
    Returns:
        (tuple, (pd.DataFrame, pd.DataFrame)): combined count and meta
            dataframes (count, meta).
    """
    count = combine_count_data([x['counts'] for x in datasets])
    meta = combine_meta_data([x['meta'] for x in datasets])
    return (count, meta)


def combine_count_data(count_matrices):
    """
    Merge count matrices on gene names.

    Merge count matrices on gene names (rows). Performs an outer join and fills
    missing values with zeros, as these genes would have been removed from a
    specific count matrix if no samples had non-zero mappings.

    Args:
        count_matrices (list, pd.DataFrame): count matrices from each
            independent scPipe run.
    Returns:
        (pd.DataFrame): merged count data from all runs. 
    """
    for i in range(1, len(count_matrices)):
        count_matrices[0] = count_matrices[0].merge(count_matrices[i],
                                                    how='outer',
                                                    left_index=True,
                                                    right_index=True)
    
    return count_matrices[0].fillna(0, inplace=True)


def combine_meta_data(meta_matrices):
    """Combine meta matrices"""

    return pd.concat(meta_matrices, axis=0, sort=True)

if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False

    if snakemake_exists:
        # only merge if flagged
        msg = "No data merging performed."
        if snakemake.params['flag']:
            # save un-merged data
            rewrite_originals(snakemake.input['cmat'], snakemake.input['meta'])

            # merge count and meta data
            outdir = os.path.commonpath([snakemake.input['cmat'],
                                         snakemake.input['meta']])
            datasets = load_datasets(snakemake.params['dirs'] + [outdir])
            combined = combine_datasets(datasets)

            # write data to original input locations 
            combined[0].to_csv(snakemake.input['cmat'])
            combined[1].to_csv(snakemake.input['meta'])

            # write output file expected by snakemake
            msg = "Succesfully combined datasets.\n\n"\
                  + "Merged data can be found in the original locations: \n"\
                  + "\tcount: {}\n\tmeta: {}\n".format(snakemake.input['cmat'],
                                                       snakemake.input['meta'])\
                  + "Original, un-merged data is also housed in the same "\
                  + "original parent directory with the 'original_' prefix.\n"\
                  + "\nCombined data from {} runs with data".format(
                                                             len(datasets))\
                  + "from {} genes over {} samples.".format(combined[0].shape[0],
                                                            combined[0].shape[1])
    
        with open(snakemake.output['out_file'], 'w') as f:
            f.write(msg)

