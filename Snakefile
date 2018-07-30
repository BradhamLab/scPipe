
# retrieve config file
import os
import re

configfile: 'files/config.yaml'
DATA_DIR = config['data_dir']
OUTPUT = config['output_dir']
LOGS = config['log_dir']
SAMPLE_REGEX = config['sample_regex']
ENDS = config['end_denote']
GENOME_DIR = os.path.join(OUTPUT, 'star/genome')
GTF = '/home/dakota/SequenceData/GenomeAnnotations/lv_genome_annotations.gtf'

# function to link sample ids to their input directory
def link_sample_dirs(data_dir, sample_regex):
    id_to_dir = {}
    sample_pattern = re.compile(sample_regex)
    for sample_dir in os.listdir(data_dir):
        matched = re.search(sample_pattern, sample_dir)
        if matched is not None:
            sample_id = sample_dir[0:matched.span()[0]]
            data_loc = os.path.join(data_dir, sample_dir)
            id_to_dir[sample_id] = data_loc
    return id_to_dir

DIRNAMES = link_sample_dirs(DATA_DIR, SAMPLE_REGEX)
IDS = list(DIRNAMES.keys())

# try and set wildcards

rule all:
    input:
        expand(os.path.join(OUTPUT, 'qc/{sample}/{sample}_{end}_qc.fastq.gz'),
               sample=IDS, end=ENDS),
        os.path.join(OUTPUT, 'fastp_summary', 'report.html'),
        os.path.join(OUTPUT, 'fastp_summary', 'read_summary.csv')


# combine lanes for each read direction
rule fastq_combine:
    input:
        lambda wildcards: DIRNAMES[wildcards.sample]
    output:
        # temporary because we'll aligned to filtered data
        temp(os.path.join(OUTPUT, 'fastq/{sample}/{sample}_{end}.fastq.gz'))
    shell:
        'cat {input}/{wildcards.sample}*_{wildcards.end}_*.fastq.gz >> {output}'

# AfterQC with fastp
rule fastp_qc:
    input:
        r1=os.path.join(OUTPUT, 'fastq/{sample}/{sample}_R1.fastq.gz'),
        r2=os.path.join(OUTPUT, 'fastq/{sample}/{sample}_R2.fastq.gz')
    log:
        os.path.join(LOGS, 'fastp/{sample}.log')
    output:
        r1=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R1_qc.fastq.gz'),
        r2=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R2_qc.fastq.gz'),
        html=os.path.join(OUTPUT, 'qc/{sample}/fastp.html'),
        json=os.path.join(OUTPUT, 'qc/{sample}/fastp.json')
    shell:
        '(fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} '
        '-h {output.html} -j {output.json}) 2> {log}'

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

# Align with star
#rule star_align:
#    input:
#        r1=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R1_qc.fastq.gz'),
#        r2=os.path.join(OUTPUT, 'qc/{sample}/{sample}_R2_qc.fastq.gz'),
#        genome=GENOME_DIR,
#        gtf=GTF
#    output:
#        genome=protected(direction(input.genome)),

#    shell:
#        'STAR --genomeDir {input.genome} --readFilesIn {input.r1} {input.r2}'

# Aggregate QC results with MultiQC - Not working
#rule run_multiqc:
#    input:
#        'output/qc'
#    output:
#        'output/multiqc/multiqc_report.html'
#    log:
#        'logs/multiqc/multiqc.log'
#    shell:
#        '(conda activate multiqc; multiqc {input} -o output/multiqc/) 2> {log}'
    
