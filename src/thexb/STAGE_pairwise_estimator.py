"""
Author: Andrew Harris

Goal: Randomly choose trees from each chromosome and calculate average coverage and p-distance
per-sample and report back. 

Need to add: Return a set of parameter values that may work for the data tested
"""
import logging
import math
import os
from pathlib import Path
from pprint import pformat
import random
import shutil
import statistics

from pyfaidx import Fasta

from thexb.UTIL_checks import check_fasta

############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/pairwise_estimator.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def count_valid_bases(seq):
    """Counts valid bases in sequence. Soft masked bases valid"""
    valid_bases = ['A', 'T', 'G', 'C']
    valid_base_count = 0
    for base in seq:
        if base not in valid_bases:
            continue
        else:
            valid_base_count += 1
            continue
    return valid_base_count


def coverage_and_median(random_windows):
    """ Calculate the average coverage per-sample """
    per_sample_cov_list = dict()
    for f in random_windows:
        with Fasta(f.as_posix(), 'fasta') as fh:
            headers = [k for k in fh.keys()]
            for sample in headers:
                seq = fh[sample][:].seq
                valid_bases = count_valid_bases(seq)
                seq_coverage = valid_bases / len(seq)
                try:
                    per_sample_cov_list[sample].append(seq_coverage)
                except KeyError:
                    per_sample_cov_list[sample] = [seq_coverage]
                continue  
    # Calculate avergae coverage per-sample
    per_sample_cov = dict()
    per_sample_median = dict()
    for s in per_sample_cov_list.keys():
        per_sample_cov[s] = round(statistics.mean(per_sample_cov_list[s]), 4)
        per_sample_median[s] = round(statistics.median(per_sample_cov_list[s]), 4)
    return per_sample_cov, per_sample_median


def p_distance(random_windows, PW_REF):
    """ Calculate the average p-distance per-sample """
    collective_pdistance = dict()
    for f in random_windows:
        with Fasta(f.as_posix(), 'fasta') as fh:
            headers = [k for k in fh.keys()]
            try:
                assert PW_REF in headers
            except AssertionError:
                logger.error(f"Provided reference, {PW_REF}, was not found in file headers. Please check reference or input file and rerun")
                exit()
            headers.pop(headers.index(PW_REF))
            ref_seq = fh[PW_REF][:].seq
            for sample in headers:
                sample_mismatch_count = 0
                sample_seq = fh[sample][:].seq
                assert len(sample_seq) == len(ref_seq)
                # Count number of incorrect pairings
                for rb, sb, in zip(ref_seq, sample_seq):
                    try:
                        assert rb == sb
                        continue
                    except AssertionError:
                        sample_mismatch_count += 1
                        continue
                sample_pdistance = sample_mismatch_count / len(sample_seq)
                try:
                    collective_pdistance[sample].append(sample_pdistance)
                except KeyError:
                    collective_pdistance[sample] = [sample_pdistance]
                continue
    avg_pdistance = dict()
    per_sample_median = dict()
    for s in collective_pdistance.keys():
        avg_pdistance[s] = round(statistics.mean(collective_pdistance[s]), 4)
        per_sample_median[s] = round(statistics.median(collective_pdistance[s]), 4)
        continue
    return avg_pdistance, per_sample_median


def return_suggested_parameters(avg_coverage, p_distance_median):
    # The suggested coverage is given as 90% of the max coverage
    max_coverage = max([c for c in avg_coverage.values()])
    suggest_cov = round(max_coverage * 0.9, 4)

    # The suggested p-distance is given as 110% of the max p-distance
    max_pdist = max([d for d in p_distance_median.values()])
    suggested_pdist = round((max_pdist * 1.1), 4)

    median_pdist = statistics.median([d for d in p_distance_median.values()])
    std_pdist = statistics.stdev([d for d in p_distance_median.values()])

    # Z-score rounded up to nearest whole number
    try:
        zscore = math.ceil((max_pdist - median_pdist) / std_pdist)
    except ZeroDivisionError:
        zscore = 1
    return suggested_pdist, zscore, suggest_cov

############################### Main Function ################################
def pairwise_estimator(filtered_indir, PW_REF, PW_EST_PERCENT_CHROM, WORKING_DIR, LOG_LEVEL):
    """Uses n-randomly chosen windows to calculate average p-distance and average coverage"""
    set_logger_level(WORKING_DIR, LOG_LEVEL)

    chrom_dirs = [d for d in filtered_indir.iterdir() if d.is_dir()]
    all_random_chrom_files = []

    # Collect random files from each chromosome
    for chrom in sorted(chrom_dirs):
        # logger.info(f"Collecting {PW_EST_PERCENT_CHROM} random files from Chromosome {chrom.stem}")
        # List windows from dir + select random windows
        
        windows = [f for f in chrom.iterdir() if check_fasta(f) & ("-DROPPED" not in f.name)]
        logger.info(f"Collecting {int(len(windows) * PW_EST_PERCENT_CHROM)} windows from Chromosome {chrom.stem}")
        random_file_nums = [random.randint(0, int(len(windows) * PW_EST_PERCENT_CHROM)) for _ in range(int(len(windows) * PW_EST_PERCENT_CHROM))]
        try:
            random_chrom_windows = [windows[i] for i in random_file_nums]
            all_random_chrom_files += random_chrom_windows
        except IndexError:
            logger.info(f"Error encountered. Possibly too few files to run. Check input parameters and rerun.")
        continue
    # Ensure there are windows collected
    try:
        assert len(all_random_chrom_files) > 0
    except AssertionError:
        logger.error("*** No files were collected. Increase percentage of chromosome sampled and rerun ***")
        return
    logger.info(f"Number of random windows: {len(all_random_chrom_files):,}\n")
    # Calculate average + median coverage per sample + log
    avg_coverage, cov_median = coverage_and_median(all_random_chrom_files)
    logger.info("Average Coverage:")
    logger.info("-----------------")
    logger.info(pformat(avg_coverage))
    logger.info("-----------------\n")
    logger.info("Median Coverage:")
    logger.info("-----------------")
    logger.info(pformat(cov_median))
    logger.info("-----------------\n")
    # Calculate average p-distance per sample + log
    avg_p_distance, pdist_median = p_distance(all_random_chrom_files, PW_REF)
    logger.info("Average p-distance:")
    logger.info("-----------------")
    logger.info(pformat(avg_p_distance))
    logger.info("-----------------\n")
    logger.info("p-distance Median:")
    logger.info("-----------------")
    logger.info(pformat(pdist_median))
    logger.info("-----------------\n")
    # Output suggested input parameters
    suggested_pdist, zscore, suggest_cov = return_suggested_parameters(avg_coverage, pdist_median)
    logger.info("Suggest input parameters:")
    logger.info(f"Note! - Suggested parameters are calculated as (maximum sampled coverage * 0.9), (maximum sampled p-distance * 1.1)")
    logger.info(f"Note! - Suggested z-score is calculated based on maximum and median (rather than mean) suggested sampled p-distance")
    logger.info("-----------------")
    logger.info(f"max_pDistance_cutoff = {suggested_pdist}")
    logger.info(f"Zscore = {zscore}")
    logger.info(f"pairwise_coverage_cutoff = {suggest_cov}")
    logger.info("-----------------")
    return


