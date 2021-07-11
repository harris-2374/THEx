import logging
import os
import shutil
import subprocess
import math
from multiprocessing import Process
from shlex import quote

import time

from ete3 import Tree
import pandas as pd

############################### Set up logger #################################
logger = logging.getLogger(__name__)

def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/iq_tree.log'):
        os.remove(WORKING_DIR / 'logs/iq_tree.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/iq_tree.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def divide_input_into_cpu_size_chunks(l, n):
    """Divides chromosomes into sets of size n, where n
    is the number of cores available to use"""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def clean_fasta_input(f):
    """Ignore hidden and .fai files. Only return valid files"""
    if f.name[0] == ".":
        return False
    elif f.suffix == ".fai":
        return False
    elif f.is_file():
        return True


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
    for chrom in chrom_dirs:
        chrom_files = [f for f in chrom.iterdir() if f.is_file()]
        for f in chrom_files:
            if "-DROPPED" in f.name:
                chrom, _, end = f.stem.strip('-DROPPED.fasta').split("_")
                treeviewer_df.at[index, 'Chromosome'] = str(chrom)
                treeviewer_df.at[index, 'Window'] = int(end)
                treeviewer_df.at[index, 'NewickTree'] = str("NoTree")
                index += 1
                continue
            else:
                chrom, start, end = f.stem.strip('.fasta').split("_")
                tree = open(f).readlines()[0]
                treeviewer_df.at[index, 'Chromosome'] = str(chrom)
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
def run_iqtree(chrom, filtered_metadata_outdir, filtered_tree_outdir, IQT_MODEL, IQT_BOOTSTRAP, IQT_CORES):
    """For each file in a given chromosome directory, run it through IQ-TREE with provided parameters"""
    metadata_suffixes = ['.bionj', '.ckp.gz', '.log', '.mldist', '.contree', '.iqtree', '.splits.nex', '.lsf', '.output', '.phy']
    chrom_files = [f for f in chrom.iterdir() if clean_fasta_input(f)]
    logger.info(f"Chromosome {chrom.stem} has {len(chrom_files)} files to run")
    skipped_files = list()
    # Iterate through each file in the chrom dir
    for cf in chrom_files:
        # Run file through IQ-TREE
        try:
            model = quote(IQT_MODEL).strip("'")
            # -nt AUTO picks best number of threads to use
            subprocess.run([f'iqtree -nt {IQT_CORES} -s {quote(cf.as_posix())} -m {quote(IQT_MODEL)} -bb {IQT_BOOTSTRAP} --quiet'], shell=True, check=True)
            continue
        except:
            skipped_files.append(cf)
            logger.error(f'IQ-Tree raised an error for {cf.name} and skipped the file\n')
            chrom_dir = filtered_tree_outdir / chrom.stem
            tree_cf = chrom_dir / f"{cf.name}-DROPPED.treefile"
            chrom_dir.mkdir(parents=True, exist_ok=True)
            with open(tree_cf, 'w') as oh:
                oh.write("WINDOW DROPPED")
                continue
    # Move IQ-TREE output to respective tree + metadata directories
    # NOTE: Not sure why, but program stalls out if I try to organize as I create files
    for iqf in chrom.iterdir():
        if iqf.suffix == '.fasta':
            continue
        elif iqf.suffix == '.treefile':
                logger.debug(f"Tree file: {iqf.name}")
                chrom_dir = filtered_tree_outdir / chrom.stem
                chrom_dir.mkdir(parents=True, exist_ok=True)
                new_filehandle = chrom_dir / iqf.name
                shutil.move(iqf, new_filehandle)
                continue
        else:
            # Check metadata suffixes
            for sufx in metadata_suffixes:
                if sufx in str(iqf):
                    logger.debug(f"Metadata file: {iqf.name}")
                    chrom_dir = filtered_metadata_outdir / chrom.stem
                    chrom_dir.mkdir(parents=True, exist_ok=True)
                    new_filehandle = chrom_dir / iqf.name
                    shutil.move(iqf, new_filehandle)
                    continue
                else:
                    continue
    return


def iq_tree(filtered_indir, filtered_metadata_outdir, filtered_tree_outdir, treeViewer_filename, 
            WORKING_DIR, IQT_MODEL, IQT_BOOTSTRAP, IQT_CORES, MULTIPROCESS, LOG_LEVEL):
    """Iterate through each trimal filetered chromosome directory and filter windows based on missingness."""
    
    st = time.time()
    
    set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level

    # Set cpu count for multiprocessing
    if type(MULTIPROCESS) == int:
        # Ensure not asking for more than available
        assert int(MULTIPROCESS) <= os.cpu_count()
        # Since each job may use more than one core, 
        # divide the number of requested cores by number of cores 
        # requested. Then round down as to not ask for one too many.
        if int(MULTIPROCESS) < IQT_CORES:
            cpu_count = 1
        else:
            cpu_count = math.floor(int(MULTIPROCESS) / IQT_CORES)
    elif MULTIPROCESS == 'all':
        cpu_count = os.cpu_count()

    # Collect Chromosome firs
    chrom_dirs = [f for f in filtered_indir.iterdir() if f.is_dir()]
    chrom_sets = list(divide_input_into_cpu_size_chunks(chrom_dirs, cpu_count))

    # Run all files through IQ-Tree
    processes = []
    for cset in chrom_sets:
        processes.clear()
        for chrom in cset:
            iqtree_args = [
                chrom,
                filtered_metadata_outdir,
                filtered_tree_outdir,
                IQT_MODEL,
                IQT_BOOTSTRAP,
                IQT_CORES,
            ]
            processes.append(Process(target=run_iqtree, args=iqtree_args))

        for process in processes:
            process.start()

        for process in processes:
            process.join()
    # Update all .treefile removing heterotachy info and secondary tree's 
    update_tree_files(filtered_tree_outdir)

    # Create initial input file for Tree Viewer
    create_TreeViewer_input(filtered_tree_outdir, treeViewer_filename)
    return

