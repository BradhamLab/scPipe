"""
Create a csv containing metadata for each sample .

Extracts data from sample ids via user provided regex pattern matching.

Author: Dakota Hawkins
Date: August 16, 2018
"""
import sys
import os

import pandas as pd

from utils import data_from_sample_name

def extract_sample_data(sample_ids, data_regex):
    """
    Extract sample metadata from sample ids using regex patterns.

    Args:
        sample_ids (list, string): list of sample ids. 
        data_regex (dict): dictionary of dictionaries relating extractable
            metadata from sample names with their associated regex patterns
            (e.g. {'Treatment': {'Chlorate': '^Chlorate'}}).

    Returns:
        (pd.DataFrame): n x k data frame where n is the number of samples and
            k is the number of extractable metadata features. 
    """
    data_series = {}
    for feature in data_regex:
        if feature != 'id':
            data_dict = data_from_sample_name(sample_ids, data_regex[feature])
            data_series[feature] = pd.Series(data_dict, index=data_dict.keys(),
                                             name=feature)
    return pd.DataFrame(data_series)

if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False

    if snakemake_exists:
        matrix = pd.read_csv(snakemake.input['csv'], index_col=0,
                             nrows=0)  # just get ids
        metadata = extract_sample_data(matrix.columns.values,
                                       snakemake.params['regex'])
        metadata['run.id'] = snakemake.params['run_id']
        metadata.to_csv(snakemake.output['csv'])
        