# system level imports
import os
import sys 

# add scripts to python path for utility functions
sys.path.append('scripts/python')
import utils

# run configuration
configfile: 'files/config.yaml'
config = utils.configure_run(config)


DATA_DIR = config['dirs']['data']
DIRNAMES = utils.link_ids_to_input(config['dirs']['data'],
                                   config['sample_regex']['id'])
IDS = list(DIRNAMES.keys())

print(IDS)
# Must align reads before creating count matrix
print("Output:\n\t{}".format("\n\t".join(utils.run_output(config))))

subworkflow read_alignment:
    workdir:
        "./subroutines/alignment"
    snakefile:
        "./subroutines/alignment/Snakefile"


rule all:
    input:
        utils.run_output(config)


rule run_pipeline:
    input:
        read_alignment(expand(os.path.join(config['dirs']['output'], 'counts',
                                           '{sample}.txt'), sample=IDS))
    output:
        os.path.join(config['dirs']['output'], 'scPipe.out')
    shell:
        'echo "Alignment complete!" > {output}'


# summarize `fastp` filtered reads
# rule summarize_fastp:
#     input:
#         os.path.join(config['dirs']['output'], 'scPipe.out')
#     params:
#         fastp=os.path.join(config['dirs']['output'], 'qc'),
#         outdir=os.path.join(config['dirs']['output'], 'fastp_summary'),
#         regex=config['sample_regex']['treatment'],
#         bad=config['thresholds']['bad'],
#         ugly=config['thresholds']['ugly']
#     output:
#         os.path.join(config['dirs']['output'], 'fastp_summary', 'report.html'),
#         os.path.join(config['dirs']['output'], 'fastp_summary',
#                      'read_summary.csv')
#     script:
#         'scripts/python/summarize_read_counts.py'


# combine counts into matrix
rule create_count_matrix:
    input:
        os.path.join(config['dirs']['output'], 'scPipe.out')
    params:
        dir=os.path.join(config['dirs']['output'], 'counts')
    output:
        count=os.path.join(config['dirs']['output'], 'matrix',
                           'count_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'matrix', 'tpm_matrix.csv')
    script:
        'scripts/python/merge_read_counts.py'

# run multiqc
rule run_multiqc:
    params:
       dir=config['dirs']['output'],
       loc=os.path.join(config['dirs']['output'], 'multiqc')
    output:
       os.path.join(config['dirs']['output'], 'multiqc', 'multiqc_report.html')
    shell:
        'conda activate multiqc; multiqc {params.dir} -o {params.loc} -m fastp '
        '-m featureCounts -m star'

# extract sample metadata from sample ids
rule create_metadata:
    input:
        csv=os.path.join(config['dirs']['output'], 'matrix', 'count_matrix.csv')
    params:
        regex=config['sample_regex'],
        run_id=config['dataset']['id']
    output:
        csv=os.path.join(config['dirs']['output'], 'metadata', 'metadata.csv')
    script:
        'scripts/python/sample_metadata.py'

# combine count matrices and metadata between datasets if flagged.
# aka some hacky bullshit where we re-write count matrices and metadata
# with combined data sources so we can pretend Snakemake is dynamic
rule combine_data:
    input:
        cmat=os.path.join(config['dirs']['output'], 'matrix',
                          'count_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'matrix', 'tpm_matrix.csv'),
        meta=os.path.join(config['dirs']['output'], 'metadata', 'metadata.csv')
    params:
        flag=config['flags']['combine_data'],
        dirs=config['dirs']['matrices']
    output:
        out_file=os.path.join(config['dirs']['output'], 'combined.out')
    script:
        'scripts/python/combine_datasets.py'

# Collapse count matrix and feature annotations down to gene level from
# transcript level.
# rule collapse_annotations:
#     input:
#         flag=os.path.join(config['dirs']['output'], 'combined.out'),
#         cmat=os.path.join(config['dirs']['output'], 'matrix',
#                           'count_matrix.csv'),
#         meta=os.path.join(config['dirs']['output'], 'metadata', 'metadata.csv'),
#         annos=config['files']['gene_annos']
#     output:
#         cmat=os.path.join(config['dirs']['output'], 'matrix',
#                           'collapsed_counts.csv'),
#         annos=os.path.join(config['dirs']['output'], 'annotations',
#                            'collapsed_annotations.csv')
#     script:
#         'scripts/python/collapse_to_genes.py'


# pre-process the count matrix performing gene/sample filtering.
rule preprocess_data:
    input:
        cmat=os.path.join(config['dirs']['output'], 'matrix',
                          'count_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'matrix', 'tpm_matrix.csv'),
        meta=os.path.join(config['dirs']['output'], 'metadata', 'metadata.csv'),
        after_combine=os.path.join(config['dirs']['output'], ('combined.out'))
    params:
        reads=config['thresholds']['bad'],
        cells=config['thresholds']['cells'],
        cpm=config['thresholds']['cpm']
    output:
        cmat=os.path.join(config['dirs']['output'], 'matrix',
                          'filtered_count_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'matrix',
                         'filtered_tpm_matrix.csv'),
        meta=os.path.join(config['dirs']['output'], 'metadata',
                          'filtered_metadata.csv')
    script:
        'scripts/python/preprocess_data.py'

# normalize raw data between samples and remove batch effects.
rule normalize_data:
    input:
        cmat=os.path.join(config['dirs']['output'], 'matrix',
                          'filtered_count_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'matrix',
                         'filtered_tpm_matrix.csv'),
        meta=os.path.join(config['dirs']['output'], 'metadata',
                          'filtered_metadata.csv')
    params:
        plot_dir=os.path.join(config['dirs']['output'], 'plots')
    output:
        cmat=os.path.join(config['dirs']['output'], 'final',
                          'normalized_log_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'final',
                         'normalized_tpm_matrix.csv'),
        meta=os.path.join(config['dirs']['output'], 'final',
                          'metadata.csv')
    script:
        'scripts/r/normalize_data.R'

rule impute_dropouts:
    input:
        cmat=os.path.join(config['dirs']['output'], 'final',
                          'normalized_log_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'final',
                         'normalized_tpm_matrix.csv')
    params:
        plot_dir=os.path.join(config['dirs']['output'], 'plots')
    output:
        cmat=os.path.join(config['dirs']['output'], 'final',
                         'imputed_log_matrix.csv'),
        tpm=os.path.join(config['dirs']['output'], 'final',
                         'imputed_tpm_matrix.csv')
    script:
        'scripts/python/impute_dropouts.py'
                
