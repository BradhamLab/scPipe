import pandas as pd

def combine_reciprocal_hits(keep_df, other_df):
    """
    """
    missed_samples = set(other_df.index.values).difference(
                     set(keep_df.index.values))
    for each in missed_samples:
        hit = other_df.loc[each, 'B_id']
        if hit not in keep_df['B_id'].values:
            new_row = [hit] + [None for i in range(keep_df.shape[1] - 1)]
            keep_df.loc[each] = new_row
    return keep_df

def combine_single_hits(keep_df, other_df):
    """
    """
    new_spus = set(other_df['subject'].unique()).difference(
                                                 keep_df['B_id'].values)
    for spu in new_spus:
        scores = other_df['bitscore'][other_df['subject'] == spu]
        row = [scores.idxmax()] + [None for i in range(keep_df.shape[1] - 1)]
        keep_df.loc[spu] = row
    return keep_df


def add_uniprot_annotations(sample_df, uniprot):
    """
    """
    gene_df = pd.DataFrame(index=uniprot.index.values,
                           columns=["UniProt.ID", "UniProt.Name"],
                           dtype=str)
    for idx in uniprot.index.values:
        prot_id, prot_name = uniprot.loc[idx, 'subject'].split('|')[1:]
        if isinstance(prot_id, str) and isinstance(prot_name, str):
            gene_df.loc[idx, 'UniProt.ID'] = prot_id
            gene_df.loc[idx, 'UniProt.Name'] = prot_name
    return pd.concat([sample_df, gene_df], axis=1, join='outer', sort=False)


def add_interpro_annotations(sample_df, interpro_file):
    """
    """
    data = {'evm': [], 'IPR.IDs': [], 'IPR.Desc': []}
    with open(interpro_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            evm = line[0]
            ipr_ids = []
            desc_ids = []
            for each in line[2:]:
                ipr, desc = each.split(';')
                ipr_ids.append(ipr.strip())
                desc_ids.append(desc.strip())
            data['evm'].append(evm)
            data['IPR.IDs'].append(';'.join(ipr_ids))
            data['IPR.Desc'].append(';'.join(desc_ids))
    ipr = pd.DataFrame(data)
    ipr.set_index('evm', inplace=True)
    return pd.concat([sample_df, ipr], axis=1, join='outer', sort=False)


def add_kegg_annotations(sample_df, kegg_file):
    """
    """
    data = {'evm': [], 'KEGG.IDs': []}
    with open(kegg_file, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            data['evm'].append(line[0])
            data['KEGG.IDs'].append(line[4])
    kegg = pd.DataFrame(data)
    kegg.set_index('evm', inplace=True)
    return pd.concat([sample_df, kegg], axis=1, join='outer', sort=False)

def add_ncbi_annotations(sample_df, ncbi):
    """
    """
    gene_df = pd.DataFrame(index=uniprot.index.values,
                           columns=["NCBI.ID"], dtype=str)
    for idx in ncbi.index.values:
        gene_df.loc[idx, 'NCBI.ID'] = ncbi.loc[idx, 'subject'].split('|')[-2]
    return pd.concat([sample_df, gene_df], axis=1, join='outer', sort=False)
    
def add_trembl_annotations(sample_df, tremble):
    gene_df = pd.DataFrame(index=uniprot.index.values,
                           columns=["TrEMBL.ID"], dtype=str)
    for idx in ncbi.index.values:
        gene_df.loc[idx, 'TrEMBL.ID'] = ncbi.loc[idx, 'subject'].split('|')[1]
    return pd.concat([sample_df, gene_df], axis=1, join='outer', sort=False)

if __name__ == "__main__":
    blast_columns = ['subject', 'perc.id', 'length', 'mismatch', 'gapopen',
                     'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    protein_models = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/ProteinModels_SPU_BestHits_peptide.txt",
                                 sep='\t', index_col=0)
    transcripts_pep = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/SPU_BestHits_peptide.txt",
                                 sep='\t', index_col=0)
    transcripts_nuc = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/SPU_BestHits.txt",
                                 sep='\t', index_col=0)
    homologues = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/best_spu_aligns.blastn",
                             sep='\t', header=None, index_col=0,
                             names=blast_columns)
    uniprot = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/Echinoderm_project/sea_urchin/5. gene_function_annotation/Lytechinus_variegatus_EVM_out_pep.SwissProt.blast",
                          sep='\t', header=None, index_col=0,
                          names=blast_columns)
    interpro_file = "/home/dakota/SequenceData/GenomeAnnotations/Echinoderm_project/sea_urchin/5. gene_function_annotation/Lytechinus_variegatus_EVM_out_pep.ipr"
    kegg_file = "/home/dakota/SequenceData/GenomeAnnotations/Echinoderm_project/sea_urchin/5. gene_function_annotation/Lytechinus_variegatus_EVM_out_pep.KEGG.blast"
    ncbi = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/Echinoderm_project/sea_urchin/5. gene_function_annotation/Lytechinus_variegatus_EVM_out_pep.nr.blast",
                       sep='\t', header=None, index_col=0, names=blast_columns)
    trembl = pd.read_csv("/home/dakota/SequenceData/GenomeAnnotations/Echinoderm_project/sea_urchin/5. gene_function_annotation/Lytechinus_variegatus_EVM_out_pep.TrEMBL.blast",
                         sep='\t', header=None, index_col=0,
                         names=blast_columns)
    annotations = combine_reciprocal_hits(pd.DataFrame(protein_models['B_id']),
                                          pd.DataFrame(transcripts_pep['B_id']))
    annotations = combine_reciprocal_hits(annotations,
                                          pd.DataFrame(transcripts_nuc))
    annotations = combine_single_hits(annotations, homologues)
    annotations.columns.values[0] = 'SPU'
    annotations = add_uniprot_annotations(annotations, uniprot)
    annotations = add_interpro_annotations(annotations, interpro_file)
    annotations = add_kegg_annotations(annotations, kegg_file)
    annotations = add_ncbi_annotations(annotations, ncbi)
    annotations = add_trembl_annotations(annotations, trembl)
    annotations.to_csv('/home/dakota/SequenceData/evm_annotations.csv')