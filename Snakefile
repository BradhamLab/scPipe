import os

configfile: 'files/config.yaml'

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


SAMPLE_REGEX = config['sample_regex']
ENDS = config['end_denote']
DATA_DIR = config['data_dir']
DIRNAMES = link_sample_dirs(DATA_DIR, SAMPLE_REGEX)
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
        "scPipe.out"

rule run_pipeline:
    input:
        quality_control((os.path.join(OUTPUT, 'fastp_summary',
                                     'read_summary.csv'),
                        expand(os.path.join(OUTPUT, 'qc', '{sample}',
                                            '{sample}_{end}_qc.fastq.gz'),
                               sample=IDS, end=ENDS))),
        read_alignment((os.path.join(OUTPUT, 'matrix', 'count_matrix.csv'),
                        expand(os.path.join(OUTPUT, 'counts', '{sample}.txt'),
                               sample=IDS)))
    output:
        'scPipe.out'
    shell:
        'echo "Run complete!" > scPipe.out'
        