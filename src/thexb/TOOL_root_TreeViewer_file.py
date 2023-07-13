import logging
import os

from ete3 import Tree
from tqdm import tqdm

from thexb.UTIL_file_reader import read_file_to_df
from thexb.STAGE_topobinner import topobinner
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
def root_TreeViewer_file(INPUT, updated_TV_dir, TV_OUTGROUP, TV_OUTGROUP_REMOVE, WORKING_DIR, LOG_LEVEL):
    """Takes in a TreeViewer input file and roots the Newick Trees with a provided outgroup"""
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    logger.debug(f"Loading {INPUT.name} into memmory")
    tv_df = read_file_to_df(INPUT)
    updated_tv_df = tv_df.copy()
    logger.debug(f"Setting new outgroup")
    log_list=[]
    with tqdm(total=len(updated_tv_df["NewickTree"])) as pbar:
        for n, t in enumerate(updated_tv_df["NewickTree"]):
            if t == "NoTree":
                pbar.update(1)
                continue
            else:
                try:
                    tree = Tree(t)
                    ancestor = tree.get_common_ancestor(TV_OUTGROUP)
                    tree.set_outgroup(ancestor)
                    if (not tree.check_monophyly(TV_OUTGROUP, target_attr="name")[0]) and (not TV_OUTGROUP_REMOVE):
                        updated_tv_df.at[n, "NewickTree"] = "NoTree"
                    else:
                        updated_tv_df.at[n, "NewickTree"] = tree.write()
                    pbar.update(1)
                    continue
                except ValueError as e:
                    log_list.append(f"{e} - {tv_df.at[n, 'Chromosome']}_{tv_df.at[n, 'Window']} - {tree.write()}")
                    continue
    # Re-bin topologies
    del updated_tv_df['TopologyID']
    # updated_tv_df['TopologyID'] = ['']*len(updated_tv_df)
    updated_tv_df.insert(3, "TopologyID", ['']*len(updated_tv_df))
    output_filename = updated_TV_dir / f"{INPUT.stem}.rooted_to_{'_'.join(TV_OUTGROUP)}.tmp{INPUT.suffix}"
    output_updated_df(updated_tv_df, INPUT, output_filename)
    topobinner(
        output_filename,
        updated_TV_dir / f"{INPUT.stem}.rooted_to_{'_'.join(TV_OUTGROUP)}{INPUT.suffix}",
        True,
        WORKING_DIR,
        LOG_LEVEL,
    )
    output_filename.unlink()
    for msg in log_list:
        logger.info(msg)
    
    return