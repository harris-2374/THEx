
import pandas as pd

def parse_treeviewer_per_chromosome(INPUT, WORKING_DIR):
    """Parse whole-genome TreeViewer file by chromosome 
       and output into new files"""
    # Read in input file
    if 'csv' in INPUT.suffix:
        df = pd.read_csv(INPUT)
    elif 'tsv' in INPUT.suffix:
        df = pd.read_csv(INPUT, sep='\t')
    elif 'xls' in INPUT.suffix:
        df = pd.read_excel(INPUT, engine='openpyxl')
    
    # Group data by chromosome
    df_grouped = df.groupby(by='Chromosome')

    for chrom, data in df_grouped:
        output_filename = WORKING_DIR / f"{chrom}_{INPUT.name}.xlsx"
        data.to_excel(output_filename, index=False)
        continue
    return