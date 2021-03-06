import logging
import os
import shutil
import subprocess
import math
from multiprocessing import Manager, freeze_support
from shlex import quote
from functools import partial

from p_tqdm import p_umap
import pandas as pd

from thexb.UTIL_checks import check_fasta

############################### Set up logger #################################
logger = logging.getLogger(__name__)

def set_logger_level(WORKING_DIR, LOG_LEVEL):
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/iq_tree.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger


############################# Custom Exceptions ###############################
class Error(Exception):
    """Base class for other exceptions"""
    pass
class InvalidCoreCount(Error):
    """Raise when config file is required but not present"""
    pass

############################## Helper Functions ###############################
def divide_input_into_cpu_size_chunks(l, n):
    """Divides chromosomes into sets of size n, where n
    is the number of cores available to use"""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def remove_heterotachy_info(l, tf):
    open_brackets = [i for i, x in enumerate(l) if x == "["]
    close_brackets = [i for i, x in enumerate(l) if x == "]"]

    final_string = f'{l[:open_brackets[0]]}'

    for ob, cb in zip(open_brackets[1:], close_brackets[:-1]):
        final_string += l[cb+1:ob]
    final_string += l[close_brackets[-1]+1:]

    with open(tf, 'w') as fh:
        fh.write(final_string)
        return


def update_tree_files(filtered_tree_outdir):
    # Edit .treefiles; remove secondary trees and heterotachy info
    for tree_chrom in filtered_tree_outdir.iterdir():
        for tf in tree_chrom.iterdir():
            if "-DROPPED" in tf.name:
                continue
            with open(tf, 'r') as tfh:
                tfh_read = tfh.readlines()
                for l in tfh_read:
                    remove_heterotachy_info(l.strip(), tf)
                    break
    return


def create_TreeViewer_input(filtered_tree_outdir, treeViewer_filename):
    treeviewer_df = pd.DataFrame(columns=['Chromosome', 'Window', 'NewickTree'])
    chrom_dirs = [c for c in filtered_tree_outdir.iterdir() if c.is_dir()]
    index = 0
    for chromosome in chrom_dirs:
        chrom_files = [f for f in chromosome.iterdir() if f.is_file()]
        for f in chrom_files:
            if "-DROPPED" in f.name:
                chromosome, _, end = f.stem.strip('-DROPPED.fasta').split("_")
                treeviewer_df.at[index, 'Chromosome'] = str(chromosome)
                treeviewer_df.at[index, 'Window'] = int(end)
                treeviewer_df.at[index, 'NewickTree'] = str("NoTree")
                index += 1
                continue
            else:
                chromosome, _, end = f.stem.strip('.fasta').split("_")
                tree = open(f).readlines()[0]
                treeviewer_df.at[index, 'Chromosome'] = str(chromosome)
                treeviewer_df.at[index, 'Window'] = int(end)
                treeviewer_df.at[index, 'NewickTree'] = str(tree)
                index += 1
                continue
    treeviewer_df.sort_values(by=['Chromosome', "Window"], inplace=True)
    # Add blank TopologyID column
    treeviewer_df['TopologyID'] = [pd.NA]*len(treeviewer_df)
    treeviewer_df.to_excel(treeViewer_filename, index=False)
    return


############################### Main Function ################################
def run_iqtree(chromosome, filtered_metadata_outdir, filtered_tree_outdir, IQT_MODEL, IQT_BOOTSTRAP, IQT_CORES, return_dict):
    """For each file in a given chromosome directory, run it through IQ-TREE with provided parameters"""
    metadata_suffixes = ['.bionj', '.ckp.gz', '.log', '.mldist', '.treefile', '.iqtree', '.splits.nex', '.lsf', '.output', '.phy']
    chrom_files = [f for f in chromosome.iterdir() if check_fasta(f)]
    dropped_files = [f for f in chromosome.iterdir() if "-DROPPED" in f.name]
    if len(dropped_files) == len(chrom_files):
        all_dropped = True
    else:
        all_dropped = False
    # logger.info(f"Chromosome {chromosome.stem} has {len(chrom_files)} files to run")
    skipped_files = list()
    dropped_total = 0
    # Iterate through each file in the chromosome dir
    for f in chrom_files:
        if all_dropped:
            skipped_files.append(f)
            chrom_dir = filtered_tree_outdir / chromosome.stem
            filestem = f.name.strip("-DROPPED.fasta")
            tree_cf = chrom_dir / f"{filestem}-DROPPED.contree"
            chrom_dir.mkdir(parents=True, exist_ok=True)
            with open(tree_cf, 'w') as oh:
                oh.write("WINDOW DROPPED")
                continue
        elif "-DROPPED" in f.name:
            dropped_total += 1
            chrom_dir = filtered_tree_outdir / chromosome.stem
            filestem = f.name.strip("-DROPPED.fasta")
            tree_cf = chrom_dir / f"{filestem}-DROPPED.contree"
            chrom_dir.mkdir(parents=True, exist_ok=True)
            with open(tree_cf, 'w') as oh:
                oh.write("WINDOW DROPPED")
                continue
        # Run file through IQ-TREE
        try:
            output_prefix = filtered_tree_outdir / chromosome.name / f.stem
            # -nt AUTO picks best number of threads to use
            subprocess.run([f'iqtree -nt {IQT_CORES} -s {quote(f.as_posix())} -m {quote(IQT_MODEL)} -bb {IQT_BOOTSTRAP} -pre {output_prefix} --quiet'], stderr=subprocess.DEVNULL, shell=True, check=True)
            continue
        except:
            skipped_files.append(f)
            chrom_dir = filtered_tree_outdir / chromosome.stem
            filestem = f.name.strip("-DROPPED.fasta")
            tree_cf = chrom_dir / f"{filestem}-DROPPED.contree"
            chrom_dir.mkdir(parents=True, exist_ok=True)
            with open(tree_cf, 'w') as oh:
                oh.write("WINDOW DROPPED")
                continue
    # Move IQ-TREE output to respective tree + metadata directories
    # NOTE: Not sure why, but program stalls out if I try to organize as I create files
    filtered_chromosome_dir = filtered_tree_outdir / chromosome.stem
    filtered_chromosome_files = [f for f in filtered_chromosome_dir.iterdir()]
    for iqf in filtered_chromosome_files:
        if iqf.suffix == '.fasta':
            continue
        elif iqf.suffix == '.contree':
                logger.debug(f"Tree file: {iqf.name}")
                chrom_dir = filtered_tree_outdir / chromosome.stem
                chrom_dir.mkdir(parents=True, exist_ok=True)
                new_filehandle = chrom_dir / iqf.name
                shutil.move(iqf, new_filehandle)
                continue
        else:
            # Check metadata suffixes
            for sufx in metadata_suffixes:
                if sufx in str(iqf):
                    logger.debug(f"Metadata file: {iqf.name}")
                    chrom_dir = filtered_metadata_outdir / chromosome.stem
                    chrom_dir.mkdir(parents=True, exist_ok=True)
                    new_filehandle = chrom_dir / iqf.name
                    shutil.move(iqf, new_filehandle)
                    continue
                else:
                    continue
    # Collect log info
    log_info = [
        f"========== {chromosome.name} ==========",
        f"Total windows in chromosome: {len(chrom_files)}",
        f"Total previously dropped windows: {dropped_total}",
        f"Total windows skipped by IQ-TREE: {len(skipped_files)}",
        f"=======================================",
    ]
    return_dict[chromosome] = log_info
    return return_dict


def iq_tree(filtered_indir, filtered_metadata_outdir, filtered_tree_outdir, treeViewer_filename, 
            WORKING_DIR, IQT_MODEL, IQT_BOOTSTRAP, IQT_CORES, MULTIPROCESS, LOG_LEVEL):
    """Iterate through each trimal filetered chromosome directory and filter windows based on missingness."""
    set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level
    freeze_support() # For Windows support
    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        # Since each job may use more than one core, 
        # divide the number of requested cores by number of cores 
        # requested. Then round down as to not ask for one too many.
        try:
            if IQT_CORES == "AUTO":
                cpu_count = MULTIPROCESS
                pass
            elif int(MULTIPROCESS) < int(IQT_CORES):
                cpu_count = 1
                pass
            else:
                cpu_count = math.floor(int(MULTIPROCESS) / int(IQT_CORES))
                pass
        except ValueError:
            logger.info((f"'{IQT_CORES}' appears to be an invalid response. Valid responses include - [AUTO, integers 1-n]"))
            exit()
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()
    else:
        logger.info((f"Invalid entry for Multiprocess value - must be an integer less than or equal to total core count"))
        return

    # Collect Chromosome firs
    chrom_dirs = sorted([f for f in filtered_indir.iterdir() if f.is_dir()])

    # Run all files through IQ-Tree
    manager = Manager()
    return_dict = manager.dict()
    p_umap(
        partial(
            run_iqtree,
            filtered_metadata_outdir=filtered_metadata_outdir,
            filtered_tree_outdir=filtered_tree_outdir,
            IQT_MODEL=IQT_MODEL,
            IQT_BOOTSTRAP=IQT_BOOTSTRAP,
            IQT_CORES=IQT_CORES,
            return_dict=return_dict,
        ),
        chrom_dirs,
        **{"num_cpus": cpu_count},
    )
    # Update all .treefile removing heterotachy info and secondary tree's 
    update_tree_files(filtered_tree_outdir)
    # Create initial input file for Tree Viewer
    create_TreeViewer_input(filtered_tree_outdir, treeViewer_filename)
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
    return

