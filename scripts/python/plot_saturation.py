from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def saturation(reads, transcripts):
    return 1 - transcripts / reads

def format_df(df, name):
    df['saturation'] = saturation(df['n.reads'], df['unique.transcripts'])
    df['percent.reads'] = df['n.reads'] / df.shape[0]
    df['name'] = name
    return df

def plot_saturation(df):
    for each in df['name'].unique():
        subset = df[df['name'] == each]
        plt.plot(subset['percent.reads'], subset['saturation'])

if __name__ == "__main__":
    try:
        snakemake
    except NameError:
        snakemake = None
    if snakemake is not None:
        dfs = []
        for each in snakemake.input['csvs']:
            df = format_df(pd.read_csv(each), each.replace('.csv', ''))
            dfs.append(df)
        combined = pd.concat(dfs, axis=0, sort=False).reset_index(drop=True)
        plot_saturation(combined)
        plt.savefig(snakemake.output['png'])
        
