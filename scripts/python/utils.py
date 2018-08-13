import os
import re
import numpy as np
import subprocess as sbp

# General helper functions
# ========================

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


# STAR helper functions
# =====================

# function to get genomeChrBinNBits parameter for STAR alignment.
def estimate_STAR_ChrBinNbits(genome_file, read_length):
    """
    Estimate the `ChrBinNBits` parameter for genome indexing in STAR

    Estimate the `ChrBinNBits` parameter for genome indexing in STAR. Value
    must be estimated due to memory constraints caused by the large number
    of scaffolds present in some genomes (i.e. the LV genome). If estimation
    is unnecessary, flag `star_est_ChrBinNbits: False` in configuration file.

    Args:
        genome_file (string): path to fasta file containing genome reference
            sequences.
        read_length (int): length of reads from RNAseq experiment.

    Return:
        (int) new value for scaling RAM consumption

    References:
    https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf (p. 7)
    https://github.com/alexdobin/STAR/issues/103
    """
    len_call = 'grep -v ">" {} | wc | awk '.format(genome_file)\
               + "'{print $3-$1}'"
    n_ref_call = 'grep "^>" {} | wc -l'.format(genome_file)

    return_values = [None, None]
    for i, call in enumerate([len_call, n_ref_call]):
        p = sbp.Popen(call, stdin=sbp.PIPE, stdout=sbp.PIPE, stderr=sbp.PIPE,
                      shell=True)
        output, err = p.communicate()
        if p.returncode == 0:
            return_values[i] = int(output.strip())
        else:
            raise OSError(err)
    estimate = max([int(np.log2(return_values[0] / return_values[1])),
                    int(np.log2(read_length))])
    return min(18, estimate)


def get_star_genome_params(config_dict):
    """
    Extract parameters for genome indexing in STAR.

    Args:
        config_dict (dictionary): configuration dictionary created by snakemake
            via configfile: {file.name}
    Returns:
        (string): string of arguments to pass STAR.
    """

    star_genome_params = config_dict['star_genome_params']
    if config_dict['star_est_ChrBinsNbits'] == True:
        nbits = estimate_STAR_ChrBinNbits(config_dict['genome_fasta'],
                                          config_dict['read_length'])
        star_genome_params += ' --genomeChrBinNbits {}'.format(nbits)
    return star_genome_params