
# retrieve config file
configfile: 'files/config.yaml'

SAMPLES = ['MK886-PMCs-3-H09_S93','MK886-PMCs-3-H12_151937799']

rule all:
    input:
        expand('output/qc/{sample}/{sample}_R1_qc.fastq.gz', sample=SAMPLES),
        expand('output/qc/{sample}/{sample}_R2_qc.fastq.gz', sample=SAMPLES),
        expand('output/qc/{sample}/{sample}_qc.html', sample=SAMPLES),
        expand('output/qc/{sample}/{sample}_qc.json', sample=SAMPLES)

rule r1_fastq_combine:
    input:
        expand('data/fastq/{sample}', sample=SAMPLES)
    output:
        expand('output/fastq/{sample}/{sample}_R1.fastq.gz', sample=SAMPLES)
    shell:
        "cat data/fastq/{wildcards.sample}/*R1*fastq.gz >> {output}"

rule r2_fastq_combine:
    input:
        expand('data/fastq/{sample}', sample=SAMPLES)
    output:
        'output/fastq/{sample}/{sample}_R2.fastq.gz'
    shell:
        "cat data/fastq/{wildcards.sample}/*R2*fastq.gz >> {output}"

# AfterQc with fastp
rule fastp_qc:
    input:
        r1=expand('output/fastq/{sample}/{sample}_R1.fastq.gz', sample=config['samples']),
        r2=expand('output/fastq/{sample}/{sample}_R2.fastq.gz', sample=config['samples'])
    output:
        r1='output/qc/{sample}/{sample}_R1_qc.fastq.gz',
        r2='output/qc/{sample}/{sample}_R2_qc.fastq.gz',
        html='output/qc/{sample}/{sample}_qc.html',
        json='output/qc/{sample}/{sample}_qc.json'
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} -j {output.json}"
