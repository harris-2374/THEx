import os
from pyfaidx import Fasta

def break_alignment_to_chroms(f):
    """Use pyfaidx to break alignment into chromosomes"""
    new_dir = f.parent / 'per_chrom/'
    os.makedirs(new_dir, 777)
    os.system(f'faidx --split-files {f} -o {new_dir}')
    return

def divide_chrom_dirs_into_chunks(l, n):
    """Breaks input into to chunks of n-size for multiprocessing"""
    for i in range(0, len(l), n):
        yield l[i:i + n]