import logging
import os

from ete3 import Tree
from tqdm import tqdm

from thexb.UTIL_file_reader import read_file_to_df
############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/root_TreeViewer_file.log'):
        os.remove(WORKING_DIR / 'logs/root_TreeViewer_file.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/root_TreeViewer_file.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger
############################## Helper Functions ###############################
def output_updated_df(updated_tv_df, INPUT, output_filename):
    """Output updated dataframe to specified file type"""
    if 'xls' in INPUT.suffix:
        updated_tv_df.to_excel(output_filename, index=False)
    elif "csv" in INPUT.suffix:
        updated_tv_df.to_csv(output_filename, index=False)
    elif 'tsv' in INPUT.suffix:
        updated_tv_df.to_csv(output_filename, sep='\t', index=False)
    return


############################### Main Function ################################
def root_TreeViewer_file(INPUT, updated_TV_dir, TV_OUTGROUP, WORKING_DIR, LOG_LEVEL):
    """Takes in a TreeViewer input file and roots the Newick Trees with a provided outgroup"""
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    logger.debug(f"Loading {INPUT.name} into memmory")
    tv_df = read_file_to_df(INPUT)
    updated_tv_df = tv_df.copy()
    logger.debug(f"Setting new outgroup")
    log_list=[]
    with tqdm(total=len(updated_tv_df["NewickTree"])) as pbar:
        for n, t in enumerate(updated_tv_df["NewickTree"]):
            # printProgressBar(n, len(updated_tv_df), prefix = 'Progress:', suffix = 'Complete', length = 50)
            if t == "NoTree":
                pbar.update(1)
                continue
            else:
                try:
                    tree = Tree(t)
                    tree.set_outgroup(TV_OUTGROUP[0])
                    updated_tv_df.at[n, "NewickTree"] = tree.write()
                    pbar.update(1)
                    continue
                except ValueError as e:
                    chrom = tv_df.at[n, "Chromosome"]
                    window = tv_df.at[n, "Window"]
                    log_list.append(f"{e} - {chrom}_{window} - {tree.write()}")
                    continue
    for msg in log_list:
        logger.info(msg)
    outgroup_str = "_".join(TV_OUTGROUP)
    output_filename = updated_TV_dir / f"{INPUT.stem}_rooted_to_{outgroup_str}{INPUT.suffix}"
    output_updated_df(updated_tv_df, INPUT, output_filename)
    return