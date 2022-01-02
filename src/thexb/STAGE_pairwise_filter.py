"""
Author: Andrew Harris + Victor Mason
Note: Code is adopted from Victor Mason's original Python2 script 
and updated to Python3
"""
import logging
import math
import os
import textwrap
from multiprocessing import Process, Manager

from tqdm import tqdm
from pyfaidx import Fasta

from thexb.UTIL_checks import check_fasta

############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/pairwise_filter.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger


############################## Helper Functions ###############################
def WriteOUT(outfile, output):
    OUT = open(outfile, 'w')
    OUT.write(output)
    OUT.close()


def FormatDictionaryOfNucleotideSeqsToFasta(d, Ns):
    o = ''
    for k in Ns:
        o += '>%s\n%s\n' % (k, "\n".join(textwrap.wrap(d[k], width=60)))
    return o


def Return1IfValueIsGreaterThanCutoffForAllInDictionary(perccov, PCcutoff, countshit):
    """Calculates percent coverage for each sample and returns False if any sample is below PCcutoff"""
    boolean = 1 # gate open
    for n in perccov.keys():
        if perccov[n] < PCcutoff:
            boolean = 0 # if one individuals percent base coverage is below PCcutoff close gate.
            if not countshit.get(n):
                countshit[n] = 1
            else:
                countshit[n] += 1
    return boolean


def CalculatePercentAACoverageForEachSequenceInDictionary(newseqs, lenseqs, PW_MISSING_CHAR):
    """Calculates percent coverage for full-length window and returns a dictionary of coverage per newseq"""
    perccov = {}
    for n in newseqs.keys():
        num_invalid_bases = 0.0
        for char in newseqs[n]:
            if char != '-' and char != PW_MISSING_CHAR and char.isalpha():
                num_invalid_bases += 1.0
        perccov[n] = float(num_invalid_bases/lenseqs)  # NOTE: Keeping as float to match other input variable data types
    return perccov


def MaskAlignedSeqsUsingStartAndEndCoordinates(seqs, start, end, PW_MISSING_CHAR, PW_WINDOW_SIZE ,lenseqs):
    """Masks regions were p-distance """
    for n in seqs.keys():
        try:
            start[n]
        except:
            pass
        else:
            for s,e in zip(start[n], end[n]):
                mask = ''
                for char in seqs[n][s-1:e]:
                    if char == '-':
                        mask += char
                    else:
                        mask += PW_MISSING_CHAR
                seqs[n] = seqs[n][:s-1] + mask + seqs[n][e:]
    return seqs


def RecordValueIfGreaterThanTwoStandardDevsAwayFromMean(median, stddev, PW_ZSCORE, PW_PDIST_CUTOFF, pdist, i, PW_WINDOW_SIZE, start, end):
    for n in pdist.keys():
        if pdist[n] > (median+(PW_ZSCORE*stddev)) and pdist[n] > PW_PDIST_CUTOFF: # value must be significantly different according to PW_ZSCORE and also greater than p-distance cutoff.
            start[n] = start.get(n, []) + [i+1]
            end[n] = end.get(n, []) + [i+PW_WINDOW_SIZE]
    return start, end


def MakeNewReducedListAndExcludeKeyValues(d, PW_EXCLUDE_LIST):
    """Removes samples in exclusion list from windows (i.e. reference sample)"""
    newdict = {}
    for key in d.keys():
        if key not in PW_EXCLUDE_LIST:
            newdict[key] = d[key]
    return newdict 


def StandardDeviation(mylist):
    """Average of the squared differences of values from mylist from their mean value. 
    minus 1 in denominator because window is sample of larger population (Bessel's correction)
    """
    mn = Mean(mylist)
    variance = sum([(e-mn)**2 for e in mylist]) / (len(mylist)-1)
    return math.sqrt(variance) # stddev = sigma = sqrt of variance


def Median(mylist):
    sorts = sorted(mylist)
    length = len(sorts)
    if length % 2 != 0:
        return ((sorts[int(length / 2)] + sorts[int(length / 2 - 1)]) / 2.0)
    else:
        return sorts[int(length / 2)]


def Mean(mylist):
    return (sum(mylist)*1.0)/len(mylist)


def CalculatePairwiseDeletionDistanceForAlignedSeqsComparedToRefSeq(seqs, PW_REF, PW_EXCLUDE_LIST, PW_MISSING_CHAR):
    pdist = {}
    # Puts single taxa in list, otherwise does nothing
    xlist = list(PW_EXCLUDE_LIST) if type(PW_EXCLUDE_LIST) == str else PW_EXCLUDE_LIST
    for n in seqs.keys():
        countgap = 0.0
        countmissing = 0.0
        countmismatch = 0.0
        denominator = 0.0
        # b is query base, r is reference base
        for b, r in zip(seqs[n], seqs[PW_REF]):
            if b == '-':
                countgap += 1.0
            elif b.upper() == PW_MISSING_CHAR:
                countmissing += 1.0
            elif b.upper() != r.upper():
                countmismatch += 1.0
                denominator += 1.0
            else:
                denominator += 1.0
        if denominator != 0:
            pdist[n] = countmismatch + countmissing / denominator
        else:
            pdist[n] = 0.0000001 # initiate very small Non-zero frequency value
    return pdist, xlist


def SplitAlignedSeqsIntoWindows(seqs, numindivs, lenseqs, PW_WINDOW_SIZE, PW_STEP, PW_MISSING_CHAR, PW_REF, PW_ZSCORE, PW_PDIST_CUTOFF):
    """Splits sequences into sub windows, calculates p-distance, removes samples from exlusion list,
    records if sub-window has a p-dist greater than X stdevs from mean, 
    then masks the sequences where p-distance is too high for a given taxa
    """
    start = {} # if window pdist is > cutoff record start and end positions of alignment for indiv
    end = {}
    num_invalid = 0

    # NOTE: Check to ensure reference sample is in file
    # Required for calculating p-distance
    try:
        assert(PW_REF in seqs.keys())
    except AssertionError:
        raise AssertionError("Reference sample not in window")
    
    for i in range(0, lenseqs, PW_STEP):
        window = {}
        for n in seqs.keys(): # cut one window for all aligned sequences in seqs.
            window[n] = seqs[n][i:(i+PW_WINDOW_SIZE)]
        pdist, xlist = CalculatePairwiseDeletionDistanceForAlignedSeqsComparedToRefSeq(window, PW_REF, PW_REF, PW_MISSING_CHAR)
        redpdist = MakeNewReducedListAndExcludeKeyValues(pdist, xlist)

        # Skip windows that have 1 or 0 samples
        if len(redpdist.values()) <= 1:
            logger.debug(f"No valid data: {num_invalid}")
            num_invalid += 1
            continue
        # mean = Mean(redpdist.values())  # NOT USED
        median = Median(redpdist.values())
        stddev = StandardDeviation(redpdist.values())
        # we use median instead of mean because less influence from outliers
        start, end = RecordValueIfGreaterThanTwoStandardDevsAwayFromMean(median, stddev, PW_ZSCORE, PW_PDIST_CUTOFF, redpdist, i, PW_WINDOW_SIZE, start, end)
    
    # FORMAT: newseqs = {"sample1": "ATCT...ACCG", "sample2": "ATCT...ACCG"}
    # NOTE: Sequences returned are back to full window length (i.e., 100kb)
    newseqs = MaskAlignedSeqsUsingStartAndEndCoordinates(seqs, start, end, PW_MISSING_CHAR, PW_WINDOW_SIZE, lenseqs)
    return newseqs


def _divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def _make_seq_dict(fh):
    seqs = dict()
    seq_len = 0
    for s in fh.keys():
        seqs[s] = fh[s][:].seq
        seq_len = len(fh[s][:].seq)
    return seqs, seq_len


def run_pw_per_chromosome(
    chrom,
    filtered_outdir,
    PW_WINDOW_SIZE,
    PW_STEP,
    PW_PDIST_CUTOFF,
    PW_REF,
    PW_MIN_SEQ_LEN,
    PW_PC_CUTOFF,
    PW_ZSCORE,
    PW_EXCLUDE_LIST,
    PW_MISSING_CHAR,
    return_dict,
):
    dropped_files = list()
    files = [f for f in chrom.iterdir() if check_fasta(f)]
    init_dropped_files = [f.stem for f in files if '-DROPPED' in f.name]
    init_valid_files = [f.stem for f in files if '-DROPPED' not in f.name]
    filtered_chrom_outdir = filtered_outdir / f'{chrom.name}'
    filtered_chrom_outdir.mkdir(parents=True, exist_ok=True)

    countshit = dict()
    countperccovcutoff = 0
    countemptyalignments = 0

    tqdm_text = "#" + f"{chrom.name}"
    with tqdm(total=len(files), desc=tqdm_text) as pbar:
        for f in files:
            # If -DROPPED in filename, write blank output file
            if "-DROPPED" in str(f.stem):
                output_fn = filtered_chrom_outdir / f.name
                with open(output_fn, 'w') as oh:
                    oh.write("")
                    continue

            with Fasta(str(f), 'fasta') as fh:
                try:
                    assert(len(fh.keys()) > 0)
                    seqs, lenseqs = _make_seq_dict(fh)
                    numindivs = len(fh.keys())
                    newseqs = SplitAlignedSeqsIntoWindows(
                        seqs,
                        numindivs,
                        lenseqs,
                        PW_WINDOW_SIZE,
                        PW_STEP,
                        PW_MISSING_CHAR,
                        PW_REF,
                        PW_ZSCORE,
                        PW_PDIST_CUTOFF,
                    )
                    perccov = CalculatePercentAACoverageForEachSequenceInDictionary(newseqs, lenseqs, PW_MISSING_CHAR)
                    # If any window has a percov < PW_PC_CUTOFF then the entire window is rejected
                    boolean = Return1IfValueIsGreaterThanCutoffForAllInDictionary(perccov, PW_PC_CUTOFF, countshit)
                    if boolean == 1:
                        output = FormatDictionaryOfNucleotideSeqsToFasta(newseqs, seqs.keys())
                        outfile = filtered_chrom_outdir / f"{f.name}"
                        WriteOUT(outfile, output)
                    elif boolean == 0:
                        outfile = filtered_chrom_outdir / f"{f.stem}-DROPPED.fasta"
                        dropped_files.append(outfile)
                        WriteOUT(outfile, "")
                        countperccovcutoff += 1

                except AssertionError:
                    logger.info(f"{f.name} has no information -- Ignoring")
                    countemptyalignments += 1
            pbar.update(1)
    num_valid_files_remaining = len([i for i in filtered_chrom_outdir.iterdir() if "-DROPPED" not in i.name])
    log_info = [
        "====================================",
        f"--- Chromosome {chrom.stem} ---",
        f"Total windows in chromosome: {len(files)}",
        f"Total valid windows at start of run: {len(init_valid_files)}",
        f"Total windows dropped from previous steps: {len(init_dropped_files)}",
        f"Total valid windows remaining after pairwise filter: {num_valid_files_remaining}",
        f"Total windows dropped by pairwise coverage cutoff: {countperccovcutoff}",
        f"Total empty alignments after pairwise filter: {countemptyalignments}",
        f"Number of times a sample caused a window to be dropped (Multiple samples can count for a single window):",
        countshit,
        f"Files dropped by pairwise filtration: {len(dropped_files)} in total",
        dropped_files,
        "====================================",
    ]
    return_dict[chrom.name] = log_info
    return return_dict


############################### Main Function ################################
def pairwise_filter(
    filtered_indir,
    filtered_outdir,
    WORKING_DIR,
    PW_WINDOW_SIZE,
    PW_STEP,
    PW_PDIST_CUTOFF,
    PW_REF,
    PW_MIN_SEQ_LEN,
    PW_PC_CUTOFF,
    PW_ZSCORE,
    PW_EXCLUDE_LIST,
    PW_MISSING_CHAR,
    MULTIPROCESS,
    LOG_LEVEL,
):
    set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()

    chrom_dirs = [f for f in filtered_indir.iterdir() if f.is_dir()]
    chrom_sets = list(_divide_chunks(chrom_dirs, cpu_count))

    proc_num = 0
    if cpu_count == 1:
        for chrom in chrom_dirs:
            run_pw_per_chromosome(
                chrom,
                filtered_outdir,
                UPDATE_FREQ,
                PW_WINDOW_SIZE,
                PW_STEP,
                PW_PDIST_CUTOFF,
                PW_REF,
                PW_MIN_SEQ_LEN,
                PW_PC_CUTOFF,
                PW_ZSCORE,
                PW_EXCLUDE_LIST,
                PW_MISSING_CHAR,
                proc_num,
            )
            proc_num += 1
    else:
        manager = Manager()
        return_dict = manager.dict()
        processes = []
        # Iterate through chroms in sets
        for cset in chrom_sets:
            processes.clear()
            for chrom in cset:
                process_args = (
                    chrom,
                    filtered_outdir,
                    PW_WINDOW_SIZE,
                    PW_STEP,
                    PW_PDIST_CUTOFF,
                    PW_REF,
                    PW_MIN_SEQ_LEN,
                    PW_PC_CUTOFF,
                    PW_ZSCORE,
                    PW_EXCLUDE_LIST,
                    PW_MISSING_CHAR,
                    return_dict,
                )
                processes.append(Process(target=run_pw_per_chromosome, args=process_args))
                proc_num += 1

            for process in processes:
                process.start()

            for process in processes:
                process.join()
        # Log output information
        for c in return_dict.keys():
            log_info = return_dict[c]
            for i in log_info:
                if type(i) == str:
                    logger.info(i)
                elif type(i) == dict:
                    for n in i.keys():
                        logger.info(f'-- {n}: {i[n]}')
                        continue
                elif type(i) == list:
                    for f in i:
                        logger.info(f"-- {f.name}")
                        continue


    return None
