import sc_utils
import scanpy.api as sc
import numpy as np
import pandas as pd

def by_name(uniprot_name):
    """
    Return the protein name for a UniProt name with a '<name>_<species>' format.
    
    Parameters
    ----------
    uniprot_name : str
        Protein name with a '<name>_<species>' format.
    
    Returns
    -------
    str
        Uniprot name for specified gene.
    """

    if isinstance(uniprot_name, str):
        return uniprot_name.split('_')[0]
    return None

def collapse_annotations(column):
    """
    Collapse gene annotations among transcripts to a single string.
    
    Parameters
    ----------
    column : pd.Series
        Annotation column from gene annotation dataframe. 
    
    Returns
    -------
    str
        String with unique annotations separated by ';'.
    """

    values = column[~column.isna()]
    if len(values) > 0:
        return ';'.join(values.astype(str).unique())
    return None

def collapse(anno_df):
    """
    Collapse an annotated dataframe by gene name.

    Collapse an annotated dataframe by gene name, such that transcripts that
    map to the same gene name (regardless of species), will have their counts
    summed and their annotations appended.
    
    Parameters
    ----------
    anno_df : sc.AnnData
        Annotated dataframe with gene names found in a 'UniProt.Name' column.
    
    Returns
    -------
    tuple, (pd.DataFrame, pd.DataFrame)
        Tuple of collapsed expression and annotation dataframes. The first
        element corresponds to the expression dataframe, while the second
        element contains the collapsed gene annotation dataframe.
    """

    # Find named genes
    names = anno_df.var.apply(lambda x: by_name(x['UniProt.Name']), axis=1)
    unique_names = names.unique()
    collapsed_expr = dict()
    collapsed_anno = dict()
    for each in unique_names:
        # collapse data by gene name
        if each is not None:
            # extract expression and annotation data for all transcripts mapped
            # to a specific gene name
            transcripts = np.where(names == each)[0]
            profile = anno_df.X[:, transcripts]
            gene_info = anno_df.var.iloc[transcripts, :].copy()
            gene_info['Transcript.ID'] = gene_info.index
            # sum counts for each cell in all transcripts
            if len(transcripts) > 1:
                profile = profile.sum(axis=1)
            # flatten array for DataFrame formatting
            elif len(transcripts) == 1:
                profile = profile.flatten()
            # Collapse annotations by ';' for each transcript annotation
            gene_info = gene_info.apply(collapse_annotations, axis=0)
            # Set name to general name, not specified by species
            gene_info.loc['UniProt.Name'] = gene_info.loc['UniProt.Name'].\
                                                            split('_')[0]
            # Add data as Series to dictionaries
            profile = pd.Series(data=profile, index=anno_df.obs.index)
            collapsed_expr[each] = profile
            collapsed_anno[each] = gene_info

    # Format collapsed data to DataFrames
    collapsed_anno = pd.DataFrame(collapsed_anno).T
    collapsed_expr = pd.DataFrame(collapsed_expr)
    # Extract non-named transcripts
    no_names = anno_df.var.index[anno_df.var['UniProt.Name'].isna()]
    # Format non-named expression data
    no_name_expr = pd.DataFrame(data=anno_df[:, no_names].X,
                                columns=no_names, index=anno_df.obs.index)
    # Format non-named annotation data
    no_name_gene = anno_df.var.loc[no_names, :].copy()
    no_name_gene['Transcript.ID'] = no_name_gene.index
    # Combine named and non-named data
    expression = pd.concat([collapsed_expr, no_name_expr], axis=1)
    genes = pd.concat([collapsed_anno, no_name_gene])

    return expression, genes


def main(in_expr, out_expr, in_anno, out_anno, in_obs):
    """
    Collapse transcript counts and annotation files to the gene level. 
    
    Parameters
    ----------
    in_expr : str
        File path to raw, uncollapsed feature count matrix as a .csv file.
    out_expr : str
        Output file name for collapsed count matrix as a .csv file.
    in_anno : str
        File path to uncollapsed gene annotations for all scaffold/transcripts
        in aligned genome. Should be in the .csv format.
    out_anno : str
        Output file name for collapse annotation file as a .csv file.
    in_obs : str
        File path to metadata for all cells in current run. Should be in the
        .csv format.
    
    Returns
    -------
    None
        Write files to specified output files. 
    """
    data = sc_utils.create_annotated_df(in_expr, in_anno, in_obs)
    expr, genes = collapse(data)
    expr.to_csv(out_expr)
    genes.to_csv(out_anno)


if __name__ == "__main__":
    snakemake_exists = True
    try:
        snakemake
    except NameError:
        snakemake_exists = False
    if snakemake_exists:
        main(snakemake.input['cmat'], snakemake.output['cmat'],
             snakemake.input['annos'], snakemake.output['annos'],
             snakemake.input['meta'])