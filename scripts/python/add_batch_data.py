"""
Add batch-level information to each cell in a count matrix.
"""
import sys
import os

import pandas as pd

from utils import assign_batch

if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False

    if snakemake_exists:
        T_matrix = pd.read_csv(snakemake.input['csv'], index_col=0).T
        batch_assignment = assign_batch(T_matrix.index.values,
                                        snakemake.params['regex'])
        batch_df = pd.DataFrame(pd.Series(batch_assignment,
                                index=batch_assignment.keys(), name='Batch'))
        final_df = pd.concat((batch_df, T_matrix), axis=1)
        final_df.to_csv(snakemake.output['csv'])
        