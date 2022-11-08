import logging
import os

from ete3 import Tree
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


def create_TreeViewer_input(filtered_tree_outdir, treeViewer_filename):
    treeviewer_df = pd.DataFrame(columns=['Chromosome', 'Window', 'NewickTree'])
    chrom_dirs = [c for c in filtered_tree_outdir.iterdir() if c.is_dir()]
    if len(chrom_dirs) == 0:
        run_type = "single-dir"
    else:
        run_type = "sub-dirs"
    idx = 0

    if run_type == "sub-dirs":
        for chrom in chrom_dirs:
            chrom_files = [f for f in chrom.iterdir() if f.is_file()]
            for f in chrom_files:
                if "-DROPPED" in f.name:
                    chrom = f.stem.strip('-DROPPED.fasta').split("_")[0]
                    end = f.stem.strip('-DROPPED.fasta').split("_")[-1]
                    treeviewer_df.at[idx, 'Chromosome'] = str(chrom)
                    treeviewer_df.at[idx, 'Window'] = int(end)
                    treeviewer_df.at[idx, 'NewickTree'] = str("NoTree")
                    idx += 1
                    continue
                else:
                    chrom = f.stem.strip('-DROPPED.fasta').split("_")[0]
                    end = f.stem.strip('-DROPPED.fasta').split("_")[-1]
                    raw_tree = open(f).readlines()[0].strip()
                    trimmed_tree = remove_heterotachy_info(raw_tree)
                    try:
                        tree = Tree(trimmed_tree)
                    except:
                        print(tree, f)
                    treeviewer_df.at[idx, 'Chromosome'] = str(chrom)
                    treeviewer_df.at[idx, 'Window'] = int(end)
                    treeviewer_df.at[idx, 'NewickTree'] = str(tree.write())
                    idx += 1
                    continue
    elif run_type == "single-file":
        chrom_files = [f for f in filtered_tree_outdir.iterdir() if f.is_file()]
        for f in chrom_files:
            if "-DROPPED" in f.name:
                chrom = f.stem.strip('-DROPPED.fasta').split("_")[0]
                end = f.stem.strip('-DROPPED.fasta').split("_")[-1]
                treeviewer_df.at[idx, 'Chromosome'] = str(chrom)
                treeviewer_df.at[idx, 'Window'] = int(end)
                treeviewer_df.at[idx, 'NewickTree'] = str("NoTree")
                idx += 1
                continue
            else:
                chrom = f.stem.strip('-DROPPED.fasta').split("_")[0]
                end = f.stem.strip('-DROPPED.fasta').split("_")[-1]
                raw_tree = open(f).readlines()[0].strip()
                trimmed_tree = remove_heterotachy_info(raw_tree)
                try:
                    tree = Tree(trimmed_tree)
                except:
                    print(tree, f)
                treeviewer_df.at[idx, 'Chromosome'] = str(chrom)
                treeviewer_df.at[idx, 'Window'] = int(end)
                treeviewer_df.at[idx, 'NewickTree'] = str(tree.write())
                idx += 1
                continue
    treeviewer_df.sort_values(by=['Chromosome', "Window"], inplace=True)
    treeviewer_df.to_excel(treeViewer_filename, index=False)
    logger.info(f"TreeViewer file written to: {treeViewer_filename}")
    return


############################### Main Function ################################
def iq_tree_external(EXTERNAL_PATH, treeViewer_filename, WORKING_DIR, LOG_LEVEL):
    """Iterate through each trimal filetered chromosome directory and filter windows based on missingness."""
    set_logger_level(WORKING_DIR, LOG_LEVEL)  # Setup log file level
    create_TreeViewer_input(EXTERNAL_PATH, treeViewer_filename)
    return

