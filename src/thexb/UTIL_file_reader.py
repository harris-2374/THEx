import pandas as pd
import os

def check_input_columns(cols):
    expected_cols = {
        'Chromosome': str,
        'Window': int,
        'NewickTree': str,
        'TopologyID': str,
    }
    final_cols = []
    # Double check standard column names, change if wrong
    for key, value, col in zip(expected_cols.keys(), expected_cols.values(), cols[:4]):
        try:
            assert type(col) == value
            assert col == key
            final_cols.append(col)
        except AssertionError:
            final_cols.append(key)
    # Add additional data column names
    if len(cols[4:]) == 0:
        final_cols.append('None')
    else:
        for c in cols[4:]:
            final_cols.append(c)
    return final_cols


def read_file_to_df(file):
    """
    Load in file depending on file type, return pandas DataFrame in json format
    """
    file_stats = os.stat(file)
    # print(f'File Size in MegaBytes is {file_stats.st_size / (1024 * 1024)}')
    # Identify file type and open accordingly
    if file.suffix == '.csv':
        open_file = pd.read_csv(file, sep=',')
        cols = check_input_columns(open_file.columns.to_list())
        if "None" in cols:
            open_file['None'] = [0]*len(open_file)
        open_file.columns = cols
        open_file.sort_values(by=[cols[0], cols[2]], inplace=True) # sort by chromosome and topology
        return open_file
    elif file.suffix == '.tsv':
        open_file = pd.read_csv(file, sep='\t')
        cols = check_input_columns(open_file.columns.to_list())
        if "None" in cols:
            open_file['None'] = [0]*len(open_file)
        open_file.sort_values(by=[cols[0], cols[2]], inplace=True) # sort by chromosome and topology
        return open_file
    elif file.suffix == '.xlsx':
        open_file = pd.read_excel(file, engine='openpyxl')
        cols = check_input_columns(open_file.columns.to_list())
        if "None" in cols:
            open_file['None'] = [0]*len(open_file)
        open_file.columns = cols
        open_file.sort_values(by=[cols[0], cols[1]], inplace=True)
        open_file.reset_index(drop=True, inplace=True)
        return open_file
    else:
        return None