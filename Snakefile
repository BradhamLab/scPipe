
# retrieve config file
configfile: 'files/config.yaml'

ENDS = ['R1', 'R2']

rule all:
    input:
        expand('output/qc/{sample}/{sample}_{end}_qc.fastq.gz', sample=config['samples'], end=ENDS)

# combine lanes for each read direction
rule fastq_combine:
    input:
        'data/fastq/{sample}'
    output:
        # temporary because we'll aligned to filtered data
        temp('output/fastq/{sample}/{sample}_{end}.fastq.gz')
    shell:
        'cat {input}/{wildcards.sample}*_{wildcards.end}_*.fastq.gz > {output}'

# AfterQC with fastp
rule fastp_qc:
    input:
        r1='output/fastq/{sample}/{sample}_R1.fastq.gz',
        r2='output/fastq/{sample}/{sample}_R2.fastq.gz'
    log:
        'logs/fastp/{sample}.log'
    output:
        r1='output/qc/{sample}/{sample}_R1_qc.fastq.gz',
        r2='output/qc/{sample}/{sample}_R2_qc.fastq.gz',
        html='output/qc/{sample}/fastp.html',
        json='output/qc/{sample}/fastp.json'
    shell:
        '(fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json}) 2> {log}'

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
    
