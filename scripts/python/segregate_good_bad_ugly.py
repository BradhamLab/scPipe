"""
Move qc-filtered fastq files into "good", "bad", and "ugly" subdirectories
depending on the number of sequenced reads.


@author: Dakota Hawkins
@date: July 31, 2018
"""

# system imports
import sys
import os
from distutils import dir_util

# csv reading
import pandas as pd

def main(summary_csv, qc_dir, output_dir):
    summary_df = pd.read_csv(summary_csv, index_col=0)

    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)

    quality_dirs = {}
    for x in ['good', 'bad', 'ugly']:
        quality_dirs[x] = os.path.join(output_dir, x)
        if not os.path.exists(quality_dirs[x]):
            os.makedirs(quality_dirs[x])

    for idx in summary_df.index:
        current_dir = os.path.join(qc_dir, idx)
        new_dir = os.path.join(
                     quality_dirs[summary_df.loc[idx, 'quality'].lower()], idx)
        if os.path.exists(current_dir):
            dir_util.copy_tree(current_dir, new_dir)
        else:
            print('{} does not exist.'.format(current_dir))

if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False
    
    if snakemake_exists:
        main(snakemake.input['csv'], snakemake.params['qc_dir'],
             snakemake.params['outdir'])
