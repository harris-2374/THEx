"""
Author: Andrew Harris
Python 3.8
"""
import logging
import os
from pathlib import Path
import pandas as pd

############################### Set up logger #################################
logger = logging.getLogger(__name__)

def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/phybin.log'):
        os.remove(WORKING_DIR / 'logs/phybin.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/phybin.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################## Helper Functions ###############################
def organize_phybin_data(PHYBIN_EXT_PATH):
    def get_binFiles(PHYBIN_EXT_PATH):
        bFiles = []
        for f in PHYBIN_EXT_PATH.iterdir():
            if f.name[0] == '.':
                continue
            elif "WARNINGS" in f.name:
                numFiles = 0
                with open(f, 'r') as fh:
                    for l in fh:
                        if ".tree" in l:
                            numFiles += 1
                        else:
                            continue
                logger.warning(f"PhyBin WARNINGS file found, ensure all data has been run properly. {numFiles} will not be added to TreeViewer input file.")
                continue
            elif '.txt' == f.suffix:
                bFiles.append(f)
            else:
                continue
        return bFiles

    def get_tree_file_info(bins):
        treeFileInfo = []
        with open(bins, 'r') as fh:
            for l in fh:
                treePath = Path(l.strip())
                treeInfo = treePath.stem
                treeFileInfo.append(treeInfo)
        return treeFileInfo

    def get_chromosome_and_windows(tree):
        """Extract chromosome and start + end windows"""
        if "." in tree:
            tree = tree.split('.')
            chrom, start, end = tree.split("_")
        else:
            chrom, start, end = tree.split("_")
        return chrom, start, end


    binDict = {f"{b.stem}": get_tree_file_info(b) for b in get_binFiles(PHYBIN_EXT_PATH)}
    binDF = pd.DataFrame(columns=['Chr', 'Window', 'Bin'])
    binDf_index = 0
    for k in binDict.keys():
        for tree in binDict[k]:
            chrom, _, end = get_chromosome_and_windows(tree)
            binDF.at[binDf_index, 'Chr'] = chrom
            binDF.at[binDf_index, 'Window'] = int(end)
            binDF.at[binDf_index, 'Bin'] = k
            binDf_index += 1
    return binDF


def merge_phybin_data(phyBinDF, tv_excel):
    tv_excel['TopologyID'] = ['NoData']*len(tv_excel)
    for row in phyBinDF.itertuples(index=False):
        chrom, window, binID = row
        tv_row = tv_excel[(tv_excel['Chromosome'] == chrom) & (tv_excel['Window'] == window)]
        try:
            tv_index = list(tv_row.index)[0]
        except IndexError:
            print(chrom)
            print(tv_row)
            exit()
        tv_excel.at[tv_index, "TopologyID"] = binID
    return tv_excel


############################### Main Function ################################
def phybin(PHYBIN_EXT_PATH, TREEVIEWER_FN, WORKING_DIR, LOG_LEVEL):
    set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level
    # Step 1: Organize data into pandas DataFrame
    phyBinDF = organize_phybin_data(PHYBIN_EXT_PATH)
    # Step 2: Load in Tree Viewer excel file
    tv_excel = pd.read_excel(TREEVIEWER_FN, engine='openpyxl')
    # Step 3: Add Bin data to TopologyID
    mergedDF = merge_phybin_data(phyBinDF, tv_excel)
    # Step 4: Output updated Tree Viewer file
    outFileName = WORKING_DIR / f"{TREEVIEWER_FN.stem}_FINAL.xlsx"
    mergedDF.to_excel(outFileName, index=False)
    return