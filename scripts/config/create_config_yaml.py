"""
Create a configuration .yaml file for single-cell analysis pipeline.

Bradham Lab, July 20, 2018

Author: Dakota Hawkins
"""

import os
import re
import yaml
import sys


def extract_fastq_info(fastq_file):
    """
    Extract sample information from fastq filenames.

    Args:
        fastq_file (string): name of .fastq file. Assumes the file name follows
            a <sample_id>_L<#*>_R<#><...>.fastq<...> format.
    Returns:
        (dict): dictionary with keys 'sample', 'lane', 'read', file, to access
            sample id, lane number, read end, and file name.
    """
    # Find lane number
    lane_match = re.compile("_L[0-9]*_")
    lane_search = lane_match.search(fastq_file)
    lane = lane_search.group(0)[1:-1]

    # Find which paired-end sample aligns to 
    read_match = re.compile("_R[0-9]_")
    read = read_match.search(fastq_file).group(0)[1:-1]

    # Find sample id
    sample = fastq_file[0:lane_search.start(0)]

    return {'sample': sample, 'lane': lane, 'read': read, 'file': fastq_file}



def get_fastq_dict(fastq_dir):
    """
    Get sample dictionaries with fastq information.

    Args:
        fastq_dir (str): path to parent directory holding fastq
            files/directories.
    Returns:
        (dict): dictionary with sample names as keys and list of associated
            fastq files.
    """
    samples = dict()
    for path, dirs, files in os.walk(fastq_dir):
        fastq_files = [extract_fastq_info(x) for x in files if 'fastq' in x]
        if len(fastq_files) > 0:
            for i, each in enumerate(fastq_files):

                file_loc = os.path.join(path, each['file'])
                # new sample, create entry
                if i == 0:
                    samples[each['sample']] = [{'file': file_loc,
                                                'lane': each['lane'],
                                                'read': each['read']}]

                # append additional files to sample list
                else:
                    samples[each['sample']].append({'file': file_loc,
                                                    'lane': each['lane'],
                                                    'read': each['read']})
    return(samples)


def create_yaml(fastqc, outdir):
    """
    Create a yaml configuration file.
    
    Yaml file is used for quality control of scRNAseq reads using fastqc during
    the Snakemake pipeline.

    Args:
        fastqc (str): top-level directory for fastq data.
        outdir (str): location to write output yaml.
    Returns:
        None
    """

    if isinstance(fastqc, str) and os.path.exists(fastqc):
        sample_dict = get_fastq_dict(fastqc)
        if isinstance(outdir, str) and os.path.exists(outdir):
            outfile = os.path.join(outdir, 'config.yaml')
            with open(outfile, 'w') as f:
                yaml.dump(sample_dict, f)
        else:
            raise IOError("{} is not a valid directory".format(outdir))
    else:
        raise IOError("{} is not a valid directory.".format(fastqc))


def usage():
    outstr = '''Create an .yaml configuration file from a top-level fastq
directory and a desired output location.\n\

    Args:\n\
        fastqc (str): top-level directory for fastq data.\n\
        outdir (str): location to write output yaml.\n\

    Execution:\n\
        python create_config_yaml.py <fastqc> <outdir>\n\n'''
    print(outstr)  


if __name__ == '__main__':
    args = sys.argv
    if len(args) == 3:
        fastqc = args[1]
        outdir = args[2]
        create_yaml(fastqc, outdir)
    else:
        print("\n\nError: No fastq directory provided.\n\n")
        usage()
        sys.exit(1)

