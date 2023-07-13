import logging
import os
from multiprocessing import Pool, Value
from time import time

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from pyfaidx import Fasta
import pandas as pd
import numpy as np

from thexb.UTIL_checks import check_fasta
################################ Important Info ################################
"""
Input:
    - Single file or a directory containing multiple files.
    - Window size to calculate p-distance in
    - Threshold of missing data to drop window calculation (default: 0.75)

    File name format: ChromosomeName.fasta
    Input directory Structure:
        WholeGenomeInSingleDirectory/
            chr1.fasta
            chr2.fasta
            chr3.fasta
            ...
            ...

Info:
Single file input will return a single file output file while multi-file returns
output for each file as well as a cumulative file to put into p-Distance Tracer.

Do not need to provide .fai file, pyfaidx will create one if cannot be found.

Functionality:
    - Calculate p-distance in windows
    - Return nan for window where a sample has more than (provided threshold) missing data (i.e., >0.75)
"""
############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/pdistance_calculator.log'):
        os.remove(WORKING_DIR / 'logs/pdistance_calculator.log')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/pdistance_calculator.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def generate_windows(seq_len, WINDOW_SIZE_INT):
    """
    Generate non-overlapping sliding windows
    """
    windows = [
        (s, e) for s, e in zip(
            range(1, seq_len, WINDOW_SIZE_INT),
            range(WINDOW_SIZE_INT, seq_len+WINDOW_SIZE_INT, WINDOW_SIZE_INT)
        )
    ]
    # Change last window end position to length of sequence - required for pyfaidx
    windows = windows[:-1] + [(windows[-1][0], seq_len)]
    return windows


def make_init_df(chromosome, queries, windows, REFERENCE):
    """
    Generate initial dataframe with all data except p-distance values
    """
    # +1 is for reference
    chromosome_list = [chromosome]*(len(windows)*(len(queries)+1))
    start_positions_list = [w[0] for w in windows]*(len(queries)+1)
    end_positions_list = [w[1] for w in windows]*(len(queries)+1)
    windows_list = [w[0] for w in windows]*(len(queries)+1)
    sample_list = []

    for q in queries:
        sample_list = sample_list + [q]*len(windows)

    sample_list = sample_list + [REFERENCE]*len(windows)

    return pd.DataFrame({
        "Chromosome": chromosome_list,
        "Start": start_positions_list,
        "End": end_positions_list,
        "Window": windows_list,
        "Sample": sample_list,
        "Value": [pd.NA]*(len(windows)*(len(queries)+1)),
    })


def pairwise_pi(seq1, seq2, PDIST_MISSING_CHAR, PDIST_IGNORE_N, PDIST_THRESHOLD):
    """
    Written by Jonas Lescroart - 10/26/2022
    Edited by Andrew Harris - 10/26/2022 - 4/10/2023

    Returns a single per-site pi value for two DNA sequences of equal length.
    Sites with missing data are optionally ignored in the final calculation.
    """
    variable_sites = 0
    invariable_sites = 0
    missing_sites = 0

    for i,j in zip(seq1.upper(), seq2.upper()):
        if PDIST_MISSING_CHAR in (i, j):
            missing_sites += 1
            continue
        elif i == j:
            invariable_sites += 1
            continue
        elif i != j:
            variable_sites += 1
            continue
        else:
            raise ValueError(f"Invalid base {i} or {j}")

    if variable_sites == invariable_sites == 0:
        return pd.NA
    elif (missing_sites/(len(seq1))) >= PDIST_THRESHOLD:
        return pd.NA
    elif not PDIST_IGNORE_N:
        return (variable_sites + missing_sites)/(len(seq1))
    else:
        return variable_sites/(invariable_sites + variable_sites)


def process_file(f, WINDOW_SIZE_INT, PDIST_MISSING_CHAR, PDIST_THRESHOLD, REFERENCE, PDIST_IGNORE_N, PDIST_REF_SUFFIX):
    """
    Load fasta file and calculate p-distance for file. Return resulting dataframe.   
    """
    chromosome = str(f.stem).replace(".fasta", "").replace(".fa", "").replace(".fna", "").replace(".fas", "")
    # Load each chromosome file
    with Fasta(f) as alignment:
        queries = [i for i in alignment.keys() if i != REFERENCE]
        logger.debug(f"{f.name} alignment loaded, starting p-distance calculation")
        # Generate windows
        windows = generate_windows(len(alignment[REFERENCE][:].seq), WINDOW_SIZE_INT)
        # Log file information
        logger.info("========================")
        logger.info(f"File: {f.name}")
        logger.info(f"Number of windows: {len(windows)}")
        logger.info(f"Samples to test: {queries}")
        # Make pandas df to save results
        df = make_init_df(chromosome, queries, windows, REFERENCE)
        # Calculate p-distance - version 1 (testing)
        for n, row in enumerate(df.to_dict('records')):
            s1 = time()
            df.at[n, "Value"] = pairwise_pi(
                alignment[REFERENCE][row['Start']:row['End']].seq,
                alignment[row['Sample']][row['Start']:row['End']].seq,
                PDIST_MISSING_CHAR,
                PDIST_IGNORE_N,
                PDIST_THRESHOLD,
            )
            logger.debug(f"{n:,}/{len(df):,} windows complete for {chromosome} :: Time:{time()- s1:.2} seconds :: DF Memory: {df.memory_usage(deep=True).sum():,}") if n%10 == 0 else None
            continue
        if PDIST_REF_SUFFIX:
            df['Sample'] = df['Sample'].apply(lambda x: f'{x}_reference' if x == REFERENCE else x)
        df.drop(columns=['Start', 'End'], inplace=True)
    logger.debug(f"-- Completed {f.name} --")
    return df

############################### Main Function ################################
def pdistance_calculator(
    INPUT,
    pdistance_output_dir,
    PDIST_THRESHOLD,
    PDIST_FILENAME,
    PDIST_MISSING_CHAR,
    REFERENCE,
    WORKING_DIR,
    WINDOW_SIZE_INT,
    PDIST_IGNORE_N,
    PDIST_REF_SUFFIX,
    MULTIPROCESS,
    LOG_LEVEL,
):
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()
    # Collect input files
    if INPUT.is_file():
        files = [INPUT]
        pass
    elif INPUT.is_dir():
        files = [f for f in INPUT.iterdir() if check_fasta(f)]
    # Create the pool
    process_pool = Pool(processes=cpu_count)
    # Start processes in the pool
    dfs = process_pool.starmap(process_file, [(f, WINDOW_SIZE_INT, PDIST_MISSING_CHAR, PDIST_THRESHOLD, REFERENCE, PDIST_IGNORE_N, PDIST_REF_SUFFIX) for f in files])
    # Concat dataframes to one dataframe
    try:
        pdist_df = pd.concat(dfs, ignore_index=True)
    except ValueError:
        pass
    outfile = pdistance_output_dir / PDIST_FILENAME
    pdist_df.reset_index(drop=True, inplace=True)
    pdist_df.to_csv(outfile, sep='\t', index=False)
    return 

