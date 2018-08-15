import os
import re
import numpy as np
import subprocess as sbp

# General helper functions
# ========================

def link_ids_to_input(data_dir, sample_regex, replicate_regex=''):
    """
    Link sample ids to data files.

    Args:
        data_dir (string): parent directory containing all sample-specific
            directories/files.
        sample_regex (string): regex pattern to extract sample ids from
            directory/file names names. Pattern should match all characters
            proceding {sample.id}, such that extracting up to the match will
            result in sample id extraction (e.g. with a file name formatted
            "/foo/bar/{sample_id}{sample_regex}.input, the pattern
            {sample_regex} should be provided.)
        replicate_regex (string, optional): regex pattern to match samples with
            their associated replicate status. Use this to only return samples
            from a specific replicate run. Default behavior has no replicate
            filtering with a default value of ''.
    Returns:
        (dict, string): dictionary linking sample id to data files
        
    """
    sample_dict = {}
    sample_pattern = re.compile(sample_regex)
    replicate_pattern = re.compile(replicate_regex)
    for sample_data in os.listdir(data_dir):
        sample_match = re.search(sample_pattern, sample_data)
        replicate_match = re.search(replicate_pattern, sample_data)
        if sample_match is not None and replicate_match is not None:
            sample_id = sample_data[0:sample_match.span()[0]]
            data_loc = os.path.join(data_dir, sample_data)
            sample_dict[sample_id] = data_loc
    return sample_dict


def data_from_sample_name(sample_ids, regex_patterns):
    """
    Extract meta-data about samples from their sample id.

    Args:
        sample_ids (list, string): sample ids containing a batch-identifying
            pattern.
        regex_patterns (dict): dictionary of regex patterns, where keys are
            the meta-data value and values are the matching regex patterns (e.g.
            {'Chlorate': '^Chlorate'})
    Returns:
        (dict): dictionary linking sample ids to meta-data values. 
    """
    data_regex = {k:re.compile(v) for k, v in regex_patterns.items()}
    sample_dict = {x:None for x in sample_ids}
    for each in sample_ids:
        for k, v in data_regex.items():
            if v.search(each) is not None:
                sample_dict[each] = k
    return sample_dict
  

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

    star_genome_params = config_dict['params']['star_genome']
    if config_dict['flags']['star_est_ChrBinsNbits'] == True:
        nbits = estimate_STAR_ChrBinNbits(config_dict['files']['genome_fasta'],
                                          config_dict['read_length'])
        star_genome_params += ' --genomeChrBinNbits {}'.format(nbits)
    return star_genome_params