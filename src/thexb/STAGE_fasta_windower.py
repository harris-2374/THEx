"""
Author: Andrew Harris
Version: Python 3.8
Summary: This script takes in a directory path that contains per-chromosome
fasta files and parses them into n-bp windows. This script can work with 
scaffold level assemblies, but each scaffold is considered it's own 
entity. 
"""
import logging
from multiprocessing import Process
from pathlib import Path
import os
import time
import textwrap

from pyfaidx import Fasta

from thexb.UTIL_checks import check_fasta
from thexb.UTIL_converters import convert_window_size_to_int
from thexb.UTIL_parsers import divide_chrom_dirs_into_chunks

############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger(WORKING_DIR, LOG_LEVEL):
    # if os.path.exists(WORKING_DIR / 'logs/fasta_windower.log'):  # Remove existing log file if present
    #     os.remove(WORKING_DIR / 'logs/fasta_windower.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/fasta_windower.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################### Helper Function ################################
def parse_chromosome_into_windows(f, WORKING_DIR, WINDOW_SIZE_INT):
    """Split each chromosome file (or just file) into n-bp windows"""
    chromosome = str(f.stem)
    try:
        with Fasta(f.as_posix()) as fasta_file:
            headers = fasta_file.keys()
            windowed_outdir = WORKING_DIR / 'windowed_fastas' / f"{chromosome}"
            windowed_outdir.mkdir(parents=True, exist_ok=True)
            flag = True
            for header in list(headers):
                logger.info(f"Parsing {header} into {WINDOW_SIZE_INT} base pair windows for chromosome {chromosome}")
                try:
                    sample_seq_df = {str(header): list(fasta_file[str(header)][:].seq)}
                except KeyError:
                    break
                # Need to add 1 to capture last bits of sequence in last window that will be less than WINDOW_SIZE_INT
                number_of_out_seqs = (len(sample_seq_df[str(header)]) // WINDOW_SIZE_INT + 1)

                start_position = 1
                end_position = WINDOW_SIZE_INT

                for _ in range(number_of_out_seqs):
                    current_file_path = windowed_outdir / f"{chromosome}_{start_position}_{end_position}.fasta"
                    if flag:
                        if current_file_path.is_file():
                            os.remove(current_file_path)
                        else:
                            pass
                    with open(current_file_path, 'a') as current_file: 
                        seq = textwrap.wrap("".join(list(sample_seq_df[str(header)][start_position:end_position])), 80)
                        current_file.write(">{}\n".format(header))
                        current_file.write("{}\n".format("\n".join(seq)))
                        start_position += WINDOW_SIZE_INT
                        end_position += WINDOW_SIZE_INT
                flag = False
    except:
        logger.warning(f"Run failed with file {chromosome}")
    return


############################### Main Function ################################
def fasta_windower(MULTI_ALIGNMENT_DIR, WORKING_DIR, WINDOW_SIZE_STR, WINDOW_SIZE_INT, MULTIPROCESS, LOG_LEVEL):
    logger = set_logger(WORKING_DIR, LOG_LEVEL)

    
    # Check if MULTI_ALIGNMENT_DIR is a file or dir
    if MULTI_ALIGNMENT_DIR.is_file():
        chrom_files = [MULTI_ALIGNMENT_DIR]
    else:
        chrom_files = [f for f in MULTI_ALIGNMENT_DIR.iterdir() if check_fasta(f)]
    logger.info(f"Parsing {len(chrom_files)} files into {WINDOW_SIZE_STR} windows")
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()

    chrom_sets = list(divide_chrom_dirs_into_chunks(chrom_files, cpu_count))
    # This block of codes reads in chromosomes
    # in sets of size 'cpu_count' and runs them
    # in parallel
    processes = []
    try:
        for cset in chrom_sets:
            processes.clear()
            for chrom in cset:
                processes.append(Process(target=parse_chromosome_into_windows, args=(chrom, WORKING_DIR, WINDOW_SIZE_INT)))
            for process in processes:
                process.start()
                logger.debug(f"Started process {process}")
            for process in processes:
                process.join()
        return
    except KeyboardInterrupt:
        for p in processes:
            if p.is_alive():
                p.terminate()
                logger.debug(f"Terminating process {p.pid}")
        return



