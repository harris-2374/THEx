"""
Author: Andrew Harris
Python 3.8
"""
import logging
import os
from multiprocessing import Process

from pyfaidx import Fasta

from thexb.UTIL_checks import missing_sample_check, check_fasta
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
    dropped_files = []
    # Ensure sequences meet minimum length
    for ff in filtered_files:
        # Load Fasta file
        with Fasta(ff.as_posix(), 'fasta') as filtered_fasta:
            filtered_headers = [k for k in filtered_fasta.keys()]
            # Iterate through headers and drop if header is missing
            for k in filtered_headers:
                if len(filtered_fasta[k][:].seq) < TRIMAL_MIN_LENGTH:
                    # Remove old index fai file
                    fai_fn = f'{ff}.fai'
                    os.remove(fai_fn)
                    # Rename file with "-DROPPED" at end
                    drop_fn = ff.parents[0] / f"{ff.stem}-DROPPED.fasta"
                    os.rename(ff, drop_fn)
                    # Log the dropped file
                    dropped_files.append(ff)
                    logger.info(f'Dropped {ff.name} because sequence does not meet minimum length of {TRIMAL_MIN_LENGTH}')
                    break
    return dropped_files


def run_trimal_per_chrom(chrom, filtered_outdir, TRIMAL_THRESH, TRIMAL_MIN_LENGTH, UPDATE_FREQ):
    """Main call of Trimal function that takes a chromosome and runs each windowed file through Trimal"""
    files = [f for f in chrom.iterdir() if check_fasta(f)]
    init_file_count = len(files)
    # Make output chromosome directory
    filtered_chrom_outdir = filtered_outdir / f'{chrom.name}'
    filtered_chrom_outdir.mkdir(parents=True, exist_ok=True)
    # Run each window through Trimal
    for n, f in enumerate(files, 1):
        if n % UPDATE_FREQ == 0:  # Run tracking
            logger.info(f" --- Completed {n} of {len(files)} for {chrom.name} ---")
        file_output_name = filtered_chrom_outdir / f'{f.name}'
        os.system(f'trimal -fasta -in {f} -out {file_output_name} -gapthreshold {TRIMAL_THRESH}')
        continue
    # Collect new filtered files
    filtered_files = [f for f in filtered_chrom_outdir.iterdir() if check_fasta(f)]
    final_file_count = len(filtered_files)
    # Filter out files with sequence lengths below TRIMAL_MIN_LENGTH
    minimum_seq_length_check(filtered_files, TRIMAL_MIN_LENGTH)
    # Recollect filtered files
    filtered_files = [f for f in filtered_chrom_outdir.iterdir() if check_fasta(f)]
    # Remove taxa or window if sample if missing 
    # when compared to pre-Trimal file (files, filtered_files)
    logger.info(f"{init_file_count-final_file_count} file(s) were removed by Trimal for chromosome {chrom.name}")
    logger.info(f" --- Completed Chromosome {chrom.name} ---")
    return


############################### Main Function ################################
def trimal(unfiltered_indir, filtered_outdir, WORKING_DIR, TRIMAL_THRESH,
           TRIMAL_MIN_LENGTH, MULTIPROCESS, LOG_LEVEL, UPDATE_FREQ):
    """Entry point for Trimal that parses the input chromosome directories 
    into chunks equal to the number cores asked to be used."""
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    # Set process list and collect chromosome directories + make sets
    processes = []
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()
    # Collect chromosome dirs and put into sets
    chrom_dirs = [c for c in unfiltered_indir.iterdir() if c.is_dir()]
    chrom_sets = list(divide_chrom_dirs_into_chunks(chrom_dirs, cpu_count))
    # Iterate through each chromosome dir
    try:
        for cset in chrom_sets:
            processes.clear()
            for chrom in cset:
                processes.append(
                    Process(
                        target=run_trimal_per_chrom, 
                        args=(chrom, filtered_outdir, TRIMAL_THRESH, TRIMAL_MIN_LENGTH, UPDATE_FREQ),
                    )
                )
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
