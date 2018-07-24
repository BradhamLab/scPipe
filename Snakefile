
# retrieve config file
configfile: 'files/config.yaml'

SAMPLES = ['MK886-PMCs-3-H09_S93','MK886-PMCs-3-H12_151937799']
ENDS = ['R1', 'R2']

rule all:
    input:
        expand('output/qc/{sample}/{sample}_{end}_qc.fastq.gz', sample=SAMPLES, end=ENDS)

# combine lanes for each read direction
rule fastq_combine:
    input:
        'data/fastq/{sample}'
    output:
        'output/fastq/{sample}/{sample}_{end}.fastq.gz'
    shell:
        'cat {input}/{wildcards.sample}*_{wildcards.end}_*.fastq.gz > {output}'

# AfterQC with fastp
rule fastp_qc:
    input:
        r1='output/fastq/{sample}/{sample}_R1.fastq.gz',
        r2='output/fastq/{sample}/{sample}_R2.fastq.gz'
    output:
        r1='output/qc/{sample}/{sample}_R1_qc.fastq.gz',
        r2='output/qc/{sample}/{sample}_R2_qc.fastq.gz',
        html='output/qc/{sample}/{sample}_qc.html',
        json='output/qc/{sample}/{sample}_qc.json'
    shell:
        'fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json}'
    
