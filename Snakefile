
# retrieve config file
configfile: 'files/config.yaml'

SAMPLES = ['MK886-PMCs-3-H09_S93','MK886-PMCs-3-H12_151937799']
ENDS = ['R1', 'R2']

rule all:
    input:
        expand('output/fastq/{sample}/{sample}_{end}.fastq.gz', sample=SAMPLES, end=ENDS)

rule fastq_combine:
    input:
        'data/fastq/{sample}'
    output:
        'output/fastq/{sample}/{sample}_{end}.fastq.gz'
    shell:
        "cat {input}/{wildcards.sample}*_{wildcards.end}_*.fastq.gz > {output}"
