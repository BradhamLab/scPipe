import os
import re
import yaml


def extract_sample_info(fastq_file):
    """
    Extract sample information from fastq filenames.

    Args:
        fastq_file (string): name of .fastq file. Assumes the file name follows
            a <sample_id>_L<#*>_R<#><...>.fastq<...> format.
    Returns:
        (dict): dictionary with keys 'sample', 'lane', 'read' to access sample
            id, lane number, and read end.
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

    return {'sample': sample, 'lane': lane, 'read': read}



def get_fastq_files(fastq_dir):
    for path, dirs, files in os.walk(fastq_dir):
        fastq_files = [extract_fastq_ids(x) for x in files if 'fastq' in x]
        if len(fastq_files) > 0:
            
        