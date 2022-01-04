import logging
import os
from multiprocessing import Pool, Value

from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import pandas as pd

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
# Make window generator
def window_generator(start, stop, WINDOW_SIZE_INT):
    return (start + WINDOW_SIZE_INT), (stop + WINDOW_SIZE_INT)

def process_file(f, WINDOW_SIZE_INT, MISSING_CHAR, PDIST_THRESHOLD, PW_REF):
    """
    Load fasta file and calculate p-distance for file. Return resulting dataframe.   
    """
    # Make pandas df to save results
    pdist_df = pd.DataFrame()
    # Load each chromosome file
    alignment = AlignIO.read(f.as_posix(), 'fasta')
    logger.info(f"{f.name} alignment loaded, starting p-distance calculation")
    calculator = DistanceCalculator('identity')
    samples = [r.id for r in alignment]
    # Ensure reference sample name provide is found in file
    try:
        assert PW_REF in samples
    except AssertionError:
        raise AssertionError("Reference sample provided is not in fasta file")
    if pdist_df.empty:
        pdist_df = pd.DataFrame(columns=['Chromosome', 'Window', "Sample", "Value"])
    start = -WINDOW_SIZE_INT
    stop = 0
    while True:
        win_start, win_stop = window_generator(start, stop, WINDOW_SIZE_INT)
        window_aln = alignment[:, win_start:win_stop]
        if window_aln.get_alignment_length() == 0:
            break
        # Check missing data
        drop_samples = []
        for i in window_aln:
            sample_name = i.name
            sample_seq = i.seq
            missing_freq = sample_seq.count(MISSING_CHAR)/len(sample_seq)
            if missing_freq > PDIST_THRESHOLD:
                drop_samples.append(sample_name)
            else:
                continue
        dist_matrix = calculator.get_distance(window_aln)
        window_contents = {
            "Chromosome": [f.stem]*len(samples),
            "Window": [win_stop]*len(samples),
            "Sample":samples,
            "Value":dist_matrix[PW_REF],
        }
        window_df = pd.DataFrame(window_contents)
        window_df.loc[window_df['Sample'].isin(drop_samples), 'Value'] = pd.NA
        try:
            pdist_df = pd.concat([pdist_df, window_df])
        except ValueError:
            pass
        # Update window position
        start = win_start
        stop = win_stop
        continue
    logger.info(f"Completed {f.name}")
    return pdist_df

############################### Main Function ################################
def pdistance_calculator(INPUT, pdistance_output_dir, PDIST_THRESHOLD, PW_REF, MISSING_CHAR, WORKING_DIR, WINDOW_SIZE_INT, MULTIPROCESS, LOG_LEVEL):
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    sum_size = 0
    files = []
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()
    # Collect input files
    if INPUT.is_file():
        chromosome_files = [INPUT]
        pass
    elif INPUT.is_dir():
        chromosome_files = [f for f in INPUT.iterdir() if check_fasta(f)]

    # Get files based on pattern and their sum of size
    for file in reversed(chromosome_files):
        sum_size =sum_size + os.path.getsize(file)
        files.append(file)
    logger.info(f'files:{len(files)} - size:{sum_size:,} bytes - processes:{cpu_count}')
    # Create the pool
    process_pool = Pool(processes=cpu_count)
    # Start processes in the pool
    dfs = process_pool.starmap(process_file, [(f, WINDOW_SIZE_INT, MISSING_CHAR, PDIST_THRESHOLD, PW_REF) for f in files])
    # Concat dataframes to one dataframe
    try:
        pdist_df = pd.concat(dfs, ignore_index=True)
    except ValueError:
        pass
    outfile = pdistance_output_dir / 'Signal_Tracer_input.tsv'
    pdist_df.reset_index(drop=True, inplace=True)
    pdist_df.to_csv(outfile, sep='\t', index=False)
    return 

