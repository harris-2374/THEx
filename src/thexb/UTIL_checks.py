import os
from pyfaidx import Fasta


def check_fasta(f):
    """Ignore hidden and .fai files. Only return valid files"""
    if f.name[0] == ".":
        return False
    elif f.suffix == ".fai":
        return False
    elif f.is_file():
        return True


def missing_sample_check(init_files, filtered_files):
    """Remove file if header is missing from trimal output file"""
    for filtered_file in filtered_files:
        for init_file in init_files:
            try:
                assert init_file.name == filtered_file.name
                # Load Fasta file
                init_fasta = Fasta(init_file.as_posix(), 'fasta')
                filtered_fasta = Fasta(filtered_file.as_posix(), 'fasta')
                # List fasta headers
                init_keys = [k for k in init_fasta.keys()]
                filtered_keys = [k for k in filtered_fasta.keys()]
                # Iterate through headers and drop if header is missing
                for k in init_keys:
                    if k not in filtered_keys:
                        os.remove(filtered_file)
            except AssertionError:
                continue


