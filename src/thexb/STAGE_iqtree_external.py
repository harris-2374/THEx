import logging
import os

from ete3 import Tree
from ete3.parser.newick import NewickError
import pandas as pd

############################### Set up logger #################################
logger = logging.getLogger(__name__)

def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/iq_tree_external.log'):
        os.remove(WORKING_DIR / 'logs/iq_tree_external.log')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/iq_tree_external.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def remove_heterotachy_info(l):
    open_brackets = [i for i, x in enumerate(l) if x == "["]
    close_brackets = [i for i, x in enumerate(l) if x == "]"]
    if open_brackets and close_brackets:
        final_string = f'{l[:open_brackets[0]]}'
        for ob, cb in zip(open_brackets[1:], close_brackets[:-1]):
            final_string += l[cb+1:ob]
        final_string += l[close_brackets[-1]+1:]
        return final_string
    else:
        return l

def update_df(treeviewer_df, chrom, end, tree):
    treeviewer_df = pd.concat([
        treeviewer_df,
        pd.DataFrame({
            'Chromosome': str(chrom),
            'Window': int(end),
            'NewickTree': str(tree),
        }, index=[0]),
    ], ignore_index=True)
    return treeviewer_df


def create_TreeViewer_input(EXTERNAL_PATH, treeviewer_filename):
    treeviewer_df = pd.DataFrame(columns=['Chromosome', 'Window', 'NewickTree', 'TopologyID'])
    chrom_dirs = [c for c in EXTERNAL_PATH.iterdir() if c.is_dir()]
    run_type = "single-dir" if len(chrom_dirs) == 0 else "sub-dirs"
    if run_type == "sub-dirs":
        for chrom in chrom_dirs:
            chrom_files = [f for f in chrom.iterdir() if f.suffix in ['.contree']]
            for f in chrom_files:
                if "-DROPPED" in f.name:
                    treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                    continue
                try:
                    chrom = f.stem.strip('-DROPPED').split("_")[0]
                    end = f.stem.strip('-DROPPED').split("_")[-1]
                    raw_tree = open(f).readlines()[0].strip()
                    tree = Tree(remove_heterotachy_info(raw_tree))
                    treeviewer_df = update_df(treeviewer_df, chrom, end, tree.write())
                    continue
                except IndexError:
                    logger.info(f"File failed - Empty file: {f.name}")
                    treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                    continue
                except NewickError:
                    logger.info(f"File failed - NewickError: {f.name}")
                    treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                    continue
    elif run_type == "single-dir":
        chrom_files = [f for f in EXTERNAL_PATH.iterdir() if f.suffix in ['.contree']]
        for f in chrom_files:
            if "-DROPPED" in f.name:
                treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                continue
            try:
                chrom = f.stem.strip('-DROPPED').split("_")[0]
                end = f.stem.strip('-DROPPED').split("_")[-1]
                raw_tree = open(f).readlines()[0].strip()
                tree = Tree(remove_heterotachy_info(raw_tree))
                treeviewer_df = update_df(treeviewer_df, chrom, end, tree.write())
                continue
            except IndexError:
                logger.info(f"File failed - Empty file: {f.name}")
                treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                continue
            except NewickError:
                logger.info(f"File failed - NewickError: {f.name}")
                treeviewer_df = update_df(treeviewer_df, chrom, end, 'NoTree')
                continue
    treeviewer_df['TopologyID'] = [None]*len(treeviewer_df)
    treeviewer_df.sort_values(by=['Chromosome', "Window"], inplace=True)
    treeviewer_df.to_excel(treeviewer_filename, index=False)
    logger.info(f"TreeViewer file written to: {treeviewer_filename}")
    return


############################### Main Function ################################
def iq_tree_external(EXTERNAL_PATH, treeviewer_filename, WORKING_DIR, LOG_LEVEL):
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    create_TreeViewer_input(EXTERNAL_PATH, treeviewer_filename)
    return

