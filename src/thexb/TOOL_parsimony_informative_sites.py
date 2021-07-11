from pathlib import Path
import logging
import os
from Bio import AlignIO
import numpy as np
import pandas as pd
import plotly
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

################################ Important Info ################################
"""
1. Input is either a directory containing all windowed fasta files from entire genome
or it contains subdirectories that contain fastas for individual chromosomes

2. If a single directory is provided and it contains fasta files for all chromosomes,
it is then assumed that files follow the proper naming convention of CHR_START_STOP.fasta

3. If a single file is provided, then the program assumes you are wanting info for just that window
and will not infer any chromosome structure to the sequence.

Input:
    window file name format: chr_start_stop.fasta
    input directory input options

    WholeGenomeInSingleDirectory/
        chr1_0_100.fasta
        chr1_101_200.fasta
        chr2_0_100.fasta
        ...
        ...

    -or-

    SubDirectoryPerChromosome/
        chr1/
            chr1_0_100.fasta
            chr1_101_200.fasta
        chr2/
            chr2_0_100.fasta
        ...
        ...


Functionality
    - Select chromosome and range [start:stop] ; default = all windows
    - Return information on missing window files
    - 

"""

############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger_level(WORKING_DIR, LOG_LEVEL):
    # Remove existing log file if present
    if os.path.exists(WORKING_DIR / 'logs/parsimony_informative_site_analysis.log'):
        os.remove(WORKING_DIR / 'logs/parsimony_informative_site_analysis.log')
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/parsimony_informative_site_analysis.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger
############################## Helper Functions ###############################
def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()


def single_file_run(input_files):
    """Generate parsimony informative site summary for single window"""
    # Read through files and calculate number of parisomy informative sites 
    # per-chromosome as well as whole-genome
    assert len(input_files) == 1
    results = pd.DataFrame(columns=["Chromosome", "Start", "Stop", "NumSites", "NumCharacters"])
    for f in input_files:
        aln = AlignIO.read(open(f), "fasta")
        for i in range(len(aln[0].seq)):
            pos_data = aln[:, i]
            pos_set = list(set(pos_data))
            # sample_names = [n.id for n in aln[:, i:i]]
            print(pos_data, pos_set)
    return


def multi_chromosome_run(input_files):
    """Generate parsimony informative site summaries per-chromosome and genome wide"""
    # Read through files and calculate number of parisomy informative sites 
    # per-chromosome as well as whole-genome
    results = pd.DataFrame(columns=["Chromosome", "Start", "Stop", "NumSites", "NumCharacters"])
    for n, f in enumerate(input_files, 1):        
        printProgressBar(n, len(input_files), prefix = 'Progress:', suffix = 'Complete', length = 50)
        f = Path(f)
        chrom, start, stop = f.stem.split("_")
        num_informative_sites = 0
        num_informative_states = 0
        aln = AlignIO.read(open(f), "fasta")
        for i in range(len(aln[0].seq)):
            pos_data = aln[:, i]
            pos_set = list(set(pos_data))
            # base_counts = {base:pos_data.count(base) for base in pos_set}
            base_counts = sorted([pos_data.count(base) for base in pos_set])
            # Check frequency of missing data + gaps is not above 25%
            # of position bases... If it is > 0.25, continue 
            if (pos_data.count("N")/len(pos_data)) > 0.25:
                continue
            elif (pos_data.count("-")/len(pos_data)) > 0.25:
                continue
            # If position has informative sites, update rolling totals
            if len(pos_set) == 1:
                # Not informative, continue
                continue
            elif len(base_counts) == 2:
                if 1 in base_counts:
                    continue
                else:
                    num_informative_sites += 1
                    num_informative_states += len(pos_set)
            elif len(base_counts) > 2:
                a, b = base_counts[-2:]
                if (a == 1) or (b==1):
                    continue
                else:
                    num_informative_sites += 1
                    num_informative_states += len(pos_set)
            else:
                raise AssertionError("Something is odd and should not reach this")
        # Update results DF with info
        if num_informative_sites == 0:
            continue
        else:
            result_row = pd.DataFrame(
                {"Chromosome": str(chrom),
                "Start": int(start),
                "Stop": int(stop),
                "NumSites": int(num_informative_sites),
                "NumCharacters": int(num_informative_states)}, 
                index=[0],
            )
        if results.empty:
            results = result_row
        else:
            results = pd.concat([results, result_row], ignore_index=True)
    return results


def build_output_plots(results, parsimony_output_dir):
    results_grouped = results.groupby(by='Chromosome')
    for chrom, df in results_grouped:
        # sp = px.scatter(df, x="Start", y="NumSites")
        sp = px.line(df, x="Start", y="NumSites")
        sp.update_layout(title=f"{chrom}")
        sp.update_xaxes(title="Position")
        sp.update_yaxes(title="Number of Informative Sites")
        # Set up output filenames + write file
        outputgraphdir = parsimony_output_dir / "Chromosome_graphs"
        outputgraphdir.mkdir(parents=True, exist_ok=True)
        outputFileHTML = outputgraphdir / f"{chrom}.html"
        outputFileSVG = outputgraphdir / f"{chrom}.svg"
        sp.write_html(outputFileHTML.as_posix())
        sp.write_image(outputFileSVG.as_posix(), format="svg", height=800, width=1250)

        dc = px.density_contour(df, x="Start", y="NumSites")
        dc.update_traces(contours_coloring="fill", contours_showlabels = True)
        dc_outputFileHTML = outputgraphdir / f"{chrom}_dc_plot.html"
        dc_outputFileSVG = outputgraphdir / f"{chrom}_dc_plot.svg"
        dc.write_html(dc_outputFileHTML.as_posix())
        dc.write_image(dc_outputFileSVG.as_posix(), format="svg", height=800, width=1250)

        subplot = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02)
        subplot.add_trace(
            go.Scattergl(
                x=df['Start'],
                y=df["NumSites"],
                mode="markers",
                showlegend=False,
            ),
            row=1,
            col=1,
        )
        subplot.add_trace(
            go.Scattergl(
                x=df['Start'],
                y=df["NumSites"],
                mode="lines",
                showlegend=False,
            ),
            row=2,
            col=1,
        )
        subplot.add_trace(
            go.Histogram2dContour(
                x=df['Start'],
                y=df["NumSites"],
                colorscale = 'Jet',
                contours = dict(
                    showlabels = True,
                    labelfont = dict(
                        family = 'Raleway',
                        color = 'white'
                    )
                ),
                hoverlabel = dict(
                    bgcolor = 'white',
                    bordercolor = 'black',
                    font = dict(
                        family = 'Raleway',
                        color = 'black'
                    )
                )
            ),
            row=3,
            col=1,
        )
        subplot.update_layout(
            title=f"{chrom}"
        )
        subplot_outputFileHTML = outputgraphdir / f"{chrom}_combined_subplot.html"
        subplot_outputFileSVG = outputgraphdir / f"{chrom}_combined_subplot.svg"
        subplot.write_html(subplot_outputFileHTML.as_posix())
        subplot.write_image(subplot_outputFileSVG.as_posix(), format="svg", height=1500, width=1250)
        subplot.show()
    return


def get_stats(results):
    grouped_results = results.groupby(by="Chromosome")
    chrom_results = dict()
    for chrom, data in grouped_results:
        chrom_results[chrom] = dict(siteCount=data.NumSites.sum(), characterCount=data.NumCharacters.sum())
        chrom_results[f"{chrom}AVG"] = dict(siteCount=data.NumSites.mean(), characterCount=data.NumCharacters.mean())
        continue
    chrom_results["WholeGenome"] = dict(siteCount=results.NumSites.sum(), characterCount=results.NumCharacters.sum())
    chrom_results["WholeGenomeAVG"] = dict(siteCount=results.NumSites.mean(), characterCount=results.NumCharacters.mean())
    df = pd.DataFrame(data=chrom_results)
    return df


############################### Main Function ################################
def parsimony_informative_sites(parsimony_output_dir, INPUT, RANGE, WORKING_DIR, LOG_LEVEL):
    set_logger_level(WORKING_DIR, LOG_LEVEL)
    # Determine input type
    if INPUT.is_file():
        results_outfilename = parsimony_output_dir / f"{INPUT.stem}_parsimony_informative_site_results.tsv"
        stats_outfilename = parsimony_output_dir / f"stats.tsv"
        run_type = 'single_file'
        input_files = [INPUT]
        pass
    elif INPUT.is_dir():
        run_type = 'multi_chromosome'
        subdirs = [d for d in INPUT.iterdir() if d.is_dir()]
        if len(subdirs) > 0:  # Assumes chromsome data is organized in individual directories
            results_outfilename = parsimony_output_dir / f"parsimony_informative_site_results.tsv"
            stats_outfilename = parsimony_output_dir / f"stats.stats"
            input_files = []
            for sd in subdirs:
                sd_files = [str(f) for f in sd.iterdir() if f.suffix != ".fai"]
                input_files += sd_files
                continue
        else: # Assumed to be a single directory (chromosome)
            results_outfilename = parsimony_output_dir / f"{INPUT.name}_parsimony_informative_site_results.tsv"
            stats_outfilename = parsimony_output_dir / f"{INPUT.name}_stats.tsv"
            input_files = [str(f) for f in INPUT.iterdir() if (f.is_file()) and (f.suffix != ".fai")]

    logger.info(f"{len(input_files):,} files collected for processing")
    # Run selected run type and return count information    
    if run_type == "single_file":
        results = single_file_run(input_files)
        pass
    elif run_type == "multi_chromosome":
        # Identify file organization and output confirmation
        results = multi_chromosome_run(input_files)
    # Sort + reset index for results
    results = results.sort_values(by=["Chromosome", "Start"])
    results = results.reset_index(drop=True)
    results.to_csv(results_outfilename, sep="\t", index=False)
    # Use results provided and generate graphs + output
    build_output_plots(results, parsimony_output_dir)
    # Output summary statistics file
    print(results)
    stats = get_stats(results)
    stats.to_csv(stats_outfilename, sep="\t")    
    return






