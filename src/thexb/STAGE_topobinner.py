"""
Author: Andrew Harris
Python 3.8
"""
import logging
import os

import pandas as pd
from ete3 import Tree
from tqdm import tqdm
############################### Set up logger #################################
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    logger = logging.getLogger(__name__)
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/topobinner.log'):
        os.remove(WORKING_DIR / 'logs/topobinner.log')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/topobinner.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def remove_heterotachy_info(l):
    """Remove any information in bracketsete3 
       does not support this format of newick"""
    if ("[" not in l) and ("]" not in l):
        return l
    open_brackets = [i for i, x in enumerate(l) if x == "["]
    close_brackets = [i for i, x in enumerate(l) if x == "]"]
    final_string = f'{l[:open_brackets[0]]}'
    for ob, cb in zip(open_brackets[1:], close_brackets[:-1]):
        final_string += l[cb+1:ob]
    final_string += l[close_brackets[-1]+1:]
    return final_string

def tv_header_validation(df):
    """Return False if first four required column headers are not valid"""
    required_cols = list(df.columns[:4])
    try:
        assert required_cols == ["Chromosome", "Window", "NewickTree", "TopologyID"]
        return True
    except AssertionError:
        return False

############################### Main Function ################################
def topobinner(TREEVIEWER_FN, UPDATED_TV_FILENAME, TOPOBIN_ROOTED, WORKING_DIR, LOG_LEVEL):
    logger = set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level
    # Load in Tree Viewer excel file
    df = pd.read_excel(TREEVIEWER_FN)
    df = df.reset_index(drop=True)
    # Validate headers
    header_check = tv_header_validation(df)
    if not header_check:
        raise AssertionError("Input file headers are not valid, please ensure required headers are correct.")
    trees = df['NewickTree']
    topologies = dict()
    logger.info(f"{len(trees):,} trees to run")
    # Set root boolean value
    TOPOBIN_ROOTED == False if TOPOBIN_ROOTED == 'Y' else True
    # Bin Trees
    tqdm_text = "{}".format("topobinner").zfill(3)
    with tqdm(total=len(trees), desc=tqdm_text) as pbar:
        for n, t in enumerate(trees):
            if t == "NoTree":
                pbar.update(1)
                continue
            elif len(topologies.keys()) == 0:
                topologies[n] = {'count': 1, 'idx': [n]}
                pbar.update(1)
                continue
            else:
                new_topology = True
                for idx in topologies.keys():
                    if df.at[idx, 'NewickTree'] == "NoTree":
                        continue
                    t1 = Tree(remove_heterotachy_info(t))
                    t2 = Tree(remove_heterotachy_info(df.at[idx, 'NewickTree']))
                    comparison = t1.compare(t2, unrooted=TOPOBIN_ROOTED)
                    rf = comparison['rf']
                    if rf == 0:
                        topologies[idx]['count'] += 1
                        topologies[idx]['idx'].append(n)
                        new_topology = False
                        break
                    else:
                        continue
                if new_topology:
                    topologies[n] = {'count': 1, 'idx': [n]}
                    pbar.update(1)
                    continue
                else:
                    pbar.update(1)
                    continue
    # Sort topologies dictionary by 'count'
    topologies = {k: v for k, v in sorted(topologies.items(), key=lambda item: item[1]['count'], reverse=True)}
    # Set zfill number
    if len(topologies.keys()) < 100:
        zfillnum = 3
    elif 100 < len(topologies.keys()) < 1000:
        zfillnum = 4
    else:
        zfillnum = 5
    # Update DataFrame TopologyID column with results
    overview_df = pd.DataFrame({
        "TopologyID": [("Tree" + "{}".format(str(i)).zfill(zfillnum)) for i in range(1, len(topologies.keys())+1)],
        "Count": [topologies[i]["count"] for i in topologies.keys()],
        "Rank": [i for i in range(1, len(topologies.keys())+1)],
    })
    topoCount = 1
    for topo in topologies.keys():
        idx = topologies[topo]['idx']
        for i in idx:
            df.at[i, 'TopologyID'] = "Tree" + "{}".format(topoCount).zfill(zfillnum)
            continue
        topoCount += 1
    # Output updated Tree Viewer file
    df.to_excel(UPDATED_TV_FILENAME, index=False)
    logger.info(f"\n{overview_df}")
    return


