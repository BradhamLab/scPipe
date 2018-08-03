
# system level imports
import os
import re
import subprocess as sbp

# numerical imports
import numpy as np

# function to link sample ids to their input directory
def link_sample_dirs(data_dir, sample_regex):
    """
    Link sample ids to fastq containing folders.

    Links sampel ids to sample-specific folders containing raw .fastq.gz
    folders. 

    Args:
        data_dir (string): parent directory containing all sample-specific
            directories.
        sample_regex (string): regex pattern to extract sample ids from
            directory names.
    Returns:
        (dict, string): dictionary mapping sample ids to data directories. 
        
    """
    id_to_dir = {}
    sample_pattern = re.compile(sample_regex)
    for sample_dir in os.listdir(data_dir):
        matched = re.search(sample_pattern, sample_dir)
        if matched is not None:
            sample_id = sample_dir[0:matched.span()[0]]
            data_loc = os.path.join(data_dir, sample_dir)
            id_to_dir[sample_id] = data_loc
    return id_to_dir


# retrieve config file
configfile: 'files/config.yaml'

# set global parameter values
DATA_DIR = config['data_dir']
OUTPUT = config['output_dir']
LOGS = config['log_dir']
SAMPLE_REGEX = config['sample_regex']
ENDS = config['end_denote']
DIRNAMES = link_sample_dirs(DATA_DIR, SAMPLE_REGEX)
IDS = list(DIRNAMES.keys())

# parameters for STAR genome index creation


# ensure path to STAR genome dir exists
if not os.path.exists(os.path.dirname(config['genome_dir'])):
    os.makedirs(os.path.dirname(config['genome_dir']))

# try and set wildcards


rule all:
    input:
        expand(os.path.join(OUTPUT, 'qc/{sample}/{sample}_{end}_qc.fastq.gz'),
               sample=IDS, end=ENDS),
        directory(os.path.join(OUTPUT, 'segregated_qc', 'bad')),
        directory(os.path.join(OUTPUT, 'segregated_qc', 'good')),
        directory(os.path.join(OUTPUT, 'segregated_qc', 'ugly')),
        #directory(os.path.join(OUTPUT, 'multiqc')),
        # protected(directory(os.path.join(OUTPUT, config['genome_dir']))),
        # os.path.join(OUTPUT, 'star', '{sample}.bam')


# combine lanes for each read direction
rule fastq_combine:
    input:
        lambda wildcards: DIRNAMES[wildcards.sample]
    output:
        # temporary because we'll aligned to filtered data
        temp(os.path.join(OUTPUT, 'fastq/{sample}/{sample}_{end}.fastq.gz'))
    shell:
        'cat {input}/{wildcards.sample}*{wildcards.end}*.fastq.gz >> {output}'

# AfterQC with fastp
# TODO: clean up qc dir manually
rule fastp_qc:
    input:
        r1=os.path.join(OUTPUT, 'fastq/{sample}/{sample}_R1.fastq.gz'),
        r2=os.path.join(OUTPUT, 'fastq/{sample}/{sample}_R2.fastq.gz')
    log:
        os.path.join(LOGS, 'fastp/{sample}.log')
    params:
        p1=config['fastp_params']
    output:
        r1=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R1_qc.fastq.gz'),
        r2=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R2_qc.fastq.gz'),
        html=os.path.join(OUTPUT, 'qc/{sample}/fastp.html'),
        json=os.path.join(OUTPUT, 'qc/{sample}/fastp.json')
    shell:
        '(fastp {params.p1} -i {input.r1} -I {input.r2} -o {output.r1} -O '
        '{output.r2} -h {output.html} -j {output.json}) 2> {log}'

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




# segregate 'good', 'bad', and 'ugly' samples
# rule segregate_samples:
#     input:
#         csv=os.path.join(OUTPUT, 'fastp_summary', 'read_summary.csv'),
#         html=os.path.join(OUTPUT, 'fastp_summary', 'report.html')
#     output:
#         directory(os.path.join(OUTPUT, 'segregated_qc', 'bad')),
#         directory(os.path.join(OUTPUT, 'segregated_qc', 'good')),
#         directory(os.path.join(OUTPUT, 'segregated_qc', 'ugly'))
#     params:
#         qc_dir=os.path.join(OUTPUT, 'qc'),
#         outdir=os.path.join(OUTPUT, 'segregated_qc')
#     script:
#         'scripts/python/segregate_good_bad_ugly.py'

# Aggregate QC results with MultiQC
# rule run_multiqc:
#     params:
#         os.path.join(OUTPUT, 'segregated_qc', 'good')
#     output:
#         directory(os.path.join(OUTPUT, 'multiqc'))
#     log:
#         os.path.join(LOGS, 'multiqc', 'multiqc.log')
#     shell:
#         '(source activate multiqc; multiqc {params} -o {output}) 2> {log}'


    
