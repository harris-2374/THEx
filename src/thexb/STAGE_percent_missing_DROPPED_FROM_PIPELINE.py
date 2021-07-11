from pprint import pprint

from multiprocessing import Process
import logging
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pyfaidx import Fasta


############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    if os.path.exists(WORKING_DIR / 'logs/percent_missing.log'):  # Remove existing log file if present
        os.remove(WORKING_DIR / 'logs/percent_missing.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/percent_missing.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def _is_fasta(f):
    """Checks if file is fasta format"""
    ftypes = ['.fasta', '.fa', '.fna', '.ffn', '.faa', '.frn']
    if f.suffix in ftypes:
        return True
    else:
        logger.warning("File is not in Fasta format, MUST provide in Fasta format")
        return False


def count_valid_chars(s):
    """Count all valid A, T, G, and C's.
       Lowercase letters count as valid bases."""
    valid_bases = ['a', 't', 'c', 'g']
    valid_total = 0
    for base in s:
        if base.lower() in valid_bases:
            valid_total += 1
        else:
            continue
    return valid_total


def output_updated_fasta_seqs(f, seqs, outdir):
    """Output files that passed filtering to new directory"""
    # Setup output dir + file
    chrom_dir = f.parts[-2]
    filename = f.parts[-1]
    outdir = outdir / f"{chrom_dir}"
    outdir.mkdir(parents=True, exist_ok=True)
    out_fn = outdir / f"{filename}"

    # Write seqs to new fasta
    with open(out_fn, 'w') as fh:
        SeqIO.write(seqs, fh, "fasta")
    return


def _divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def clean_filtered_output(filtered_outdir):
    chrom_dirs = [d for d in filtered_outdir.iterdir() if d.is_dir()]
    for c in chrom_dirs:
        empty_files = [f for f in c.iterdir() if os.path.getsize(f) == 0]
        nonempty_files = [f for f in c.iterdir() if os.path.getsize(f) > 0]
        for f in nonempty_files:
            filesize = os.path.getsize(f)
            if filesize == 0:
                os.remove(f)
                print('Removed')
                continue
            else:
                continue
    return

############################### Main Function ################################
def filter_chromosome_files(WORKING_DIR, chrom, filtered_outdir, PERCENT_DROP_TYPE, PERCENT_MISSING_THRESH, UPDATE_FREQ):
    # Update console
    logger.info(f"---------- Filtering Chromosome {chrom.name} ---------- ")

    # File collection and tracking
    chrom_files = [f for f in chrom.iterdir() if (f.suffix != '.fai') and _is_fasta(f)]
    dropped_windows = []

    # Run progress tracking
    curr_file = 1

    # Iterate through each window (fasta file)
    for f in chrom_files:
        updated_seqs = []
        with Fasta(str(f), 'fasta') as ffile:
            window_drop_flag = False
            # update console progress
            if curr_file % UPDATE_FREQ == 0:
                logger.info(f"--- Completed {curr_file}/{len(chrom_files)} files for chromosome {chrom.name} ---")
            # Iterate through each sample and check sequence missing content
            for sample in ffile.keys():
                num_valid = count_valid_chars(ffile[sample][:].seq)
                seq_length = len(ffile[sample][:].seq)
                if num_valid/seq_length > PERCENT_MISSING_THRESH:
                    seqrec = SeqRecord(
                        Seq(ffile[sample][:].seq),
                        id=str(sample),
                        name=str(sample),
                        description='',  # leave empty
                    )
                    updated_seqs.append(seqrec)
                    logger.debug(f"[PASS] {f.name} - {sample} - {num_valid}/{seq_length} or {((num_valid/seq_length)*100):2f}% valid bases")
                    pass
                else:
                    if PERCENT_DROP_TYPE == 'taxa':
                        # Ignores taxa; not added to updated_seqs
                        continue
                    elif PERCENT_DROP_TYPE == 'window':
                        window_drop_flag = True
                        logger.debug(f"[DROP_WINDOW] {chrom.name} - {f.name}")
                        dropped_windows.append(str(f.name))
                        break
            # If PERCENT_DROP_TYPE=='taxa', output file with updated seqs to new file
            if PERCENT_DROP_TYPE == "taxa":
                output_updated_fasta_seqs(f, updated_seqs, filtered_outdir)
                curr_file += 1
                continue
            else:
                if window_drop_flag:
                    curr_file += 1
                    continue
                else:
                    output_updated_fasta_seqs(f, updated_seqs, filtered_outdir)
                    curr_file += 1
                    continue
    
    clean_filtered_output(filtered_outdir)
    
    # Log number of dropped windows
    if PERCENT_DROP_TYPE == 'window':
        dropped_dir = WORKING_DIR / 'dropped_files'
        dropped_dir.mkdir(parents=True, exist_ok=True)
        outfilename = dropped_dir / f'{chrom.name}_dropped_files.tsv'
        with open(outfilename, 'w') as oh:
            oh.write("\n".join(dropped_windows))
        logger.info(f"--- {len(dropped_windows)} of {len(chrom_files)} windows dropped from chromosome {chrom.name}")
        return
    else:
        return
    

def percent_missing(filtered_indir, filtered_outdir, WORKING_DIR, PERCENT_MISSING_THRESH, PERCENT_DROP_TYPE, MULTIPROCESS, LOG_LEVEL, UPDATE_FREQ):
    """Iterate through each trimal filetered chromosome directory and filter windows based on missingness."""
    logger = set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level

    # Set cpu count for multiprocessing
    if MULTIPROCESS.isnumeric():
        # Will not use more than total cpu's available - 1
        assert int(MULTIPROCESS) <= (os.cpu_count() - 1)
        cpu_count = int(MULTIPROCESS)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()

    # Collect Chromosome first
    chrom_dirs = [f for f in filtered_indir.iterdir() if f.is_dir()]

    chrom_sets = list(_divide_chunks(chrom_dirs, cpu_count))

    processes = []

    # Iterate through each chromosome dir
    try:
        for cset in chrom_sets:
            processes.clear()
            for chrom in cset:
                processes.append(Process(target=filter_chromosome_files, args=(WORKING_DIR, chrom, filtered_outdir, PERCENT_DROP_TYPE, PERCENT_MISSING_THRESH, UPDATE_FREQ)))
            for process in processes:
                process.start()
            for process in processes:
                process.join()
        return
    except KeyboardInterrupt:
        for p in processes:
            if p.is_alive():
                p.terminate()
                logger.debug(f"Terminating process {p.pid}")
        return
    return None

# if __name__ == "__main__":
#     percent_missing()