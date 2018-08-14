# system level imports
import os
import sys 

# add scripts to python path for utility functions
sys.path.append('scripts/python')
from utils import link_ids_to_input

configfile: 'files/config.yaml'

SAMPLE_REGEX = config['sample_regex']
ENDS = config['end_denote']
DATA_DIR = config['data_dir']
DIRNAMES = link_ids_to_input(DATA_DIR, SAMPLE_REGEX)
IDS = list(DIRNAMES.keys())
OUTPUT = config['output_dir']

subworkflow quality_control:
    workdir: './subroutines/qc'
    snakefile: './subroutines/qc/Snakefile'

subworkflow read_alignment:
    workdir: './subroutines/alignment'
    snakefile: './subroutines/alignment/Snakefile'

rule all:
    input:
        os.path.join(OUTPUT, 'scPipe.out'),
        os.path.join(OUTPUT, 'fastp_summary', 'report.html'),
        os.path.join(OUTPUT, 'fastp_summary', 'read_summary.csv'),
        os.path.join(OUTPUT, 'matrix', 'count_matrix.csv'),
        os.path.join(OUTPUT, 'multiqc', 'multiqc_report.html')


rule run_pipeline:
    input:
        read_alignment(expand(os.path.join(OUTPUT, 'counts', '{sample}.txt'),
                              sample=IDS))
    output:
        os.path.join(OUTPUT, 'scPipe.out')
    shell:
        'echo "Run complete!" > {output}'


# summarize `fastp` filtered reads
rule summarize_fastp:
    params:
        fastp=os.path.join(OUTPUT, 'qc'),
        outdir=os.path.join(OUTPUT, 'fastp_summary'),
        regex=config['treatment_regex'],
        bad=config['bad_threshold'],
        ugly=config['ugly_threshold']
    output:
        os.path.join(OUTPUT, 'fastp_summary', 'report.html'),
        os.path.join(OUTPUT, 'fastp_summary', 'read_summary.csv')
    script:
        'scripts/python/summarize_read_counts.py'


# combine counts into matrix
rule create_count_matrix:
    params:
        dir=os.path.join(OUTPUT, 'counts')
    output:
        csv=os.path.join(OUTPUT, 'matrix', 'count_matrix.csv')
    script:
        'scripts/python/merge_read_counts.py'

# run multiqc
rule run_multiqc:
    params:
        dir=OUTPUT,
        loc=os.path.join(OUTPUT, 'multiqc')
    output:
        os.path.join(OUTPUT, 'multiqc', 'multiqc_report.html')
    shell:
        'source activate multiqc; multiqc {params.dir} -o {params.loc} -m fastp '
        '-m featureCounts -m star'
        