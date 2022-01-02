"""
Author: Andrew Harris
Python 3.8
"""
import logging
import os
import subprocess
from multiprocessing import Process, Manager
from shlex import quote

from tqdm import tqdm
from pyfaidx import Fasta

from thexb.UTIL_checks import check_fasta
from thexb.UTIL_parsers import divide_chrom_dirs_into_chunks

############################### Set up logger #################################
logger = logging.getLogger(__name__)


def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/trimal.log'):
        os.remove(WORKING_DIR / 'logs/trimal.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/trimal.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger


############################## Helper Functions ###############################
def minimum_seq_length_check(filtered_files, TRIMAL_MIN_LENGTH):
    """Removes file from filtered output directory if output file contains sequence of lengths
    less than TRIMAL_MIN_LENGTH."""
    log_info = list()
    # Ensure sequences meet minimum length
    for f in filtered_files:
        if "-DROPPED" in f.name:
            continue
        # Load Fasta file
        with Fasta(f.as_posix(), 'fasta') as filtered_fasta:
            filtered_headers = [k for k in filtered_fasta.keys()]
            # Iterate through headers and drop if header is missing
            for k in filtered_headers:
                if len(filtered_fasta[k][:].seq) < TRIMAL_MIN_LENGTH:
                    # Remove old index fai file
                    fai_fn = f'{f}.fai'
                    os.remove(fai_fn)
                    # Rename file with "-DROPPED" at end
                    drop_fn = f.parents[0] / f"{f.stem}-DROPPED.fasta"
                    os.rename(f, drop_fn)
                    log_info.append(f'*File Dropped* | File: {f.name} | Reason: Sequence does not meet minimum length of {TRIMAL_MIN_LENGTH}')
                    break
    return log_info


def empty_seq_log(f):
    return f"*File Dropped* | File: {f} | Reason: All alignment sequence was removed by Trimal"


def write_empty_files(empty_files, filtered_chrom_outdir):
    """Write -DROPPED.fasta files for empty alignment files

    :param empty_files: list of filenames that are empty
    :type empty_files: list
    """
    for f in empty_files:
        out_file_name = filtered_chrom_outdir / f"{f.stem}-DROPPED.fasta"
        with open(out_file_name, 'w') as oh:
            oh.write("")
            continue
    return


def run_trimal_per_chrom(chrom, filtered_outdir, TRIMAL_THRESH, TRIMAL_MIN_LENGTH, return_dict):
    """Main call of Trimal function that takes a chromosome and runs each windowed file through Trimal"""
    files = [f for f in chrom.iterdir() if check_fasta(f)]
    init_file_count = len(files)
    # Make output chromosome directory
    filtered_chrom_outdir = filtered_outdir / f'{chrom.name}'
    filtered_chrom_outdir.mkdir(parents=True, exist_ok=True)
    # Run each window through Trimal
    tqdm_text = "#" + f"{chrom.name}"
    with tqdm(total=len(files), desc=tqdm_text) as pbar:
        for f in files:
            file_output_name = filtered_chrom_outdir / f'{f.name}'
            subprocess.run([f'trimal -fasta -in {quote(f.as_posix())} -out {quote(file_output_name.as_posix())} -gapthreshold {TRIMAL_THRESH}'], shell=True, check=True, stderr=subprocess.DEVNULL)
            pbar.update(1)
            continue
    # Collect new filtered files
    filtered_files = [f for f in filtered_chrom_outdir.iterdir() if check_fasta(f)]
    # Filter out files with sequence lengths below TRIMAL_MIN_LENGTH
    log_info = minimum_seq_length_check(filtered_files, TRIMAL_MIN_LENGTH)
    # Recollect filtered files
    filtered_files = [f for f in filtered_chrom_outdir.iterdir() if check_fasta(f)]
    # Identify files with no remaining sequence
    empty_files = [f for f in files if f.name not in [i.name for i in filtered_files]]
    write_empty_files(empty_files, filtered_chrom_outdir)
    # Generate log messages for empty sequences
    empty_seq_logs = [empty_seq_log(f.name) for f in empty_files]
    # Calculate remaining files
    final_valid_file_count = init_file_count - (len(log_info) + len(empty_files))
    final_fail_seq_len_file_count = len(log_info)
    final_dropped_file_count = len(empty_files)
    final_valid_file_count = init_file_count - (len(log_info) + len(empty_files))
    return_dict[chrom] = (log_info, init_file_count, final_valid_file_count, final_fail_seq_len_file_count, final_dropped_file_count, empty_seq_logs)
    return return_dict


############################### Main Function ################################
def trimal(unfiltered_indir, filtered_outdir, WORKING_DIR, TRIMAL_THRESH,
           TRIMAL_MIN_LENGTH, MULTIPROCESS, LOG_LEVEL):
    """Entry point for Trimal that parses the input chromosome directories 
    into chunks equal to the number cores asked to be used."""
    set_logger_level(WORKING_DIR, LOG_LEVEL)

    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()
    # Collect chromosome dirs and put into sets
    chrom_dirs = sorted([c for c in unfiltered_indir.iterdir() if c.is_dir()])
    chrom_sets = list(divide_chrom_dirs_into_chunks(chrom_dirs, cpu_count))
    # Iterate through each chromosome dir
    manager = Manager()
    return_dict = manager.dict()
    processes = []
    try:
        for cset in chrom_sets:
            processes.clear()
            for chrom in cset:
                processes.append(
                    Process(
                        target=run_trimal_per_chrom, 
                        args=(
                            chrom,
                            filtered_outdir,
                            TRIMAL_THRESH, 
                            TRIMAL_MIN_LENGTH,
                            return_dict,
                        ),
                    )
                )
            for process in processes:
                process.start()
                logger.debug(f"Started process {process}")
            for process in processes:
                process.join()
        
        # Log output information
        for c in return_dict.keys():
            log_info, init_file_count, final_valid_file_count, final_fail_seq_len_file_count, final_dropped_file_count, empty_seq_logs = return_dict[c]
            print()
            logger.info(f"====================== {c.name} ======================")
            logger.info("-------------------")
            logger.info(f"Inital file count: {init_file_count}")
            logger.info(f"Failed to minimum sequence length: {final_fail_seq_len_file_count}")
            logger.info(f"No sequence remaining in file: {final_dropped_file_count}")
            logger.info(f"Remaining valid files: {final_valid_file_count}")
            logger.info("-------------------")
            for i in log_info:
                logger.info(i)
            for j in empty_seq_logs:
                logger.info(j)
            logger.info(f"=====================================================")
            
        
        return
    except KeyboardInterrupt:
        for p in processes:
            if p.is_alive():
                p.terminate()
                logger.debug(f"Terminating process {p.pid}")
        return
