"""
Script to impute dropouts in normalized scRNAseq data using MAGIC

Reference:

van Dijk, D., Sharma, R., Nainys, J., Yim, K., Kathail, P., Carr, A. J., … 
Pe’er, D. (2018). Recovering Gene Interactions from Single-Cell Data Using
Data Diffusion. Cell, 174(3), 716–729.e27.
https://doi.org/10.1016/J.CELL.2018.05.061

Author: Dakota Hawkins
Date: August 24, 2018
"""

import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'  # make cluster not crash

import sys

import pandas as pd
import magic

import seaborn as sns
from matplotlib import pyplot as plt
# switch to backend for scc runs
plt.switch_backend('agg')

CLUSTER = False

if __name__ == '__main__':
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False
    if snakemake_exists:
        data_holder = {'count': {'input': snakemake.input['cmat'],
                                 'output': snakemake.output['cmat']},
                       'tpm': {'input': snakemake.input['tpm'],
                               'output': snakemake.output['tpm']}}
        
        for key in data_holder.keys():
            # read in data
            data = pd.read_csv(data_holder[key]['input'], index_col=0)
        
            # impute with magic
            magic_op = magic.MAGIC()
            imputed = magic_op.fit_transform(data.T, genes='all_genes')

            # write data
            imputed.to_csv(data_holder[key]['output'])
        
            if not CLUSTER: 
                # plot non-imputed data
                orig_heatmap = sns.clustermap(data.T, z_score=1, cmap='Blues')
                plt.savefig(os.path.join(snakemake.params['plot_dir'],
                            '{}_heatmap.png'.format(key)))
                plt.cla()

                # plot imputed data
                imputed_heatmap = sns.clustermap(imputed, z_score=1,
                                                 cmap='Blues')
                plt.savefig(os.path.join(snakemake.params['plot_dir'],
                            '{}_imputed_heatmap.png'.format(key)))
                plt.cla()

