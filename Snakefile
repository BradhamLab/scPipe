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


# Must align reads before creating count matrix
subworkflow read_alignment:
    workdir: './subroutines/alignment'
    snakefile: './subroutines/alignment/Snakefile'

rule all:
    input:
        os.path.join(config['dirs']['output'], ('combined.out'))
        # os.path.join(config['dirs']['output'], 'scPipe.out'),
        # os.path.join(config['dirs']['output'], 'fastp_summary', 'report.html'),
        # os.path.join(config['dirs']['output'], 'fastp_summary', 'read_summary.csv'),
        # os.path.join(config['dirs']['output'], 'matrix', 'count_matrix.csv'),
        # os.path.join(config['dirs']['output'], 'multiqc', 'multiqc_report.html')


rule run_pipeline:
    input:
        read_alignment(expand(os.path.join(config['dirs']['output'], 'counts',
                                           '{sample}.txt'), sample=IDS))
    output:
        os.path.join(config['dirs']['output'], 'scPipe.out')
    shell:
        'echo "Alignment complete!" > {output}'


# summarize `fastp` filtered reads
rule summarize_fastp:
    params:
        fastp=os.path.join(config['dirs']['output'], 'qc'),
        outdir=os.path.join(config['dirs']['output'], 'fastp_summary'),
        regex=config['sample_regex']['treatment'],
        bad=config['thresholds']['bad'],
        ugly=config['thresholds']['ugly']
    output:
        os.path.join(config['dirs']['output'], 'fastp_summary', 'report.html'),
        os.path.join(config['dirs']['output'], 'fastp_summary',
                     'read_summary.csv')
    script:
        'scripts/python/summarize_read_counts.py'


# combine counts into matrix
rule create_count_matrix:
    params:
        dir=os.path.join(config['dirs']['output'], 'counts')
    output:
        csv=os.path.join(config['dirs']['output'], 'matrix', 'count_matrix.csv')
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
        'source activate multiqc; multiqc {params.dir} -o {params.loc} -m fastp '
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
        meta=os.path.join(config['dirs']['output'], 'metadata', 'metadata.csv')
    params:
        flag=config['flags']['combine_data'],
        dirs=config['dirs']['matrices']
    output:
        out_file=os.path.join(config['dirs']['output'], ('combined.out'))
    script:
        'scripts/python/combine_datasets.py'



        