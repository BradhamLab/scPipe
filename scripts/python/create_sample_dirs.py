"""
Script to find all files associated with a given sample in a provided directory,
create sample specific sub-directories within the provided directory, and
finally, move sample associated files to their respecitive sub-directory.

@author: Dakota Hawkins
@date: July 31, 2018
"""

# system level imports
import os
import sys
import shutil

# regex to match sample files
import re

def relate_samples(fastq_dir):
    """
    Relate all files associated with the same sample id.

    Args:
        fastq_dir (string): path to folder containing mixed sample fastq files.
            Files names are assumed to follow a <sample_id>_<read>.fastq format.
    Returns:
        (dict, list): dictionary containing lists of all files associated with
            each provided sample. Keys are sample ids. 
    """

    files = os.listdir(fastq_dir)
    sample_ids = set(['_'.join(x.split('_')[0:-1]) for x in files\
                     if os.path.isfile(os.path.join(fastq_dir, x))])
    sample_files = {}

    for x in sample_ids:
        sample_regex = re.compile('^' + x)
        matched_files = []
        for fastq in files:
            matched = sample_regex.search(fastq)
            if matched is not None and\
            os.path.isfile(os.path.join(fastq_dir, fastq)):
                matched_files.append(os.path.join(fastq_dir, fastq))
        sample_files[x] = matched_files
    
    return sample_files


def create_sample_dirs(fastq_dir):
    """
    Create subdirectories and move all files associated with a given sample.
    """

    related_samples = relate_samples(fastq_dir)

    for key in related_samples:
        sample_dir = os.path.join(fastq_dir, key)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        for f in related_samples[key]:
            filename = os.path.basename(f)
            new_loc = os.path.join(sample_dir, filename)
            shutil.move(f, new_loc)
        
if __name__ == '__main__':
    if len(sys.argv) == 2 and isinstance(sys.argv[1], str):
        if os.path.exists(sys.argv[1]):
            create_sample_dirs(sys.argv[1])
