# Tree House Explorer

 Tree House Explorer (THEx) is a novel phylogenomic genome browser that host a growing collection of interactive dashboards. It was developed using [Plotly's](https://plotly.com/) open source framework [Dash](https://dash.plotly.com/) for building extremely powerful and customizable analytical web applications.

## THEx on Windows

Two of the main dependencies required for the THExBuilder's Tree Viewer pipeline are not supported on Windows. Due to this limitation, THEx for Windows is only deployed with the main THEx application.

## Browser Support

THEx is supported by most major browsers (Google Chrome, Safari, Firefox, etc.) except Internet Explorer.  

# Conda Installation

Steps:

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
   - Windows 10 installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)
   - Mac OS installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)  
   - Linux installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
2. Add required channels
    > $ conda config --add channels bioconda
    >
    > $ conda config --add channels conda-forge
3. Create a new conda environment
    > $ conda create -n thex_env
4. Activate conda environment
    > $ conda activate thex_env
5. Install THEx
    > $ conda install thex -c ajharris_2374  

# THEx Start-up

To launch THEx, activate your THEx conda environment and run.

```
thex
```

```
usage: thex [--host HOST] [--port PORT]

optional arguments:
  --host HOST  Host address (default: 127.0.0.1)
  --port PORT  Port number (default: 8050)
```

Once the session is running, the message ```Tree House Explorer running on http://127.0.0.1:8050/``` will appear. Click or copy/paste the address into your browser to access your session.

# THExBuilder

## Usage

THExBuilder is a built-in command line suite that contains tools and pipelines that help simplify the process of building and manipulating input files for THEx. There are currently two main pipelines, one for Tree Viewer and one for Signal Tracer. To use THExBuilder, activate your THEx conda environment and enter the command _"thexb -h"_. This command will bring up the help section, providing an overview of all the options and available tools.

```
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print version

general inputs:
  -i , --input          Input file path (Cannot be used when running Tree Viewer pipeline with configuration file)
  -o , --output         Output directory (Cannot be used when running Tree Viewer pipeline with configuration file)
  -r , --reference      Reference sample for the --pdistance, --pw_estimator, and --pw_filter commands.
  -w , --window_size    Window size [bp/kb/mb/gb] (default: 100kb)
  -c , --config         Pathway to configuration file for Tree Viewer pipeline (*recommended approach*)

Tree Viewer pipeline stages:
  --tv_all              Run all steps of the Tree Viewer pipeline sequentially
  --minifastas          Fasta windowing step
  --trimal              Run windowed fasta files through Trimal
  --pw_estimator        Additional tool that provides basic information to help choose initial settings for the pairwise filtration step
  --pw_filter           Sliding window approach to masking segments of windows where a sample's divergence is above the provided z-score threshold as well as dropping windows that contain samples with coverage below given threshold.
  --iqtree              Run windowed fasta files through IQTree
  --iqtree_external     Create base TreeViewer input file with externally run IQ-Tree results. Use (-i) to provide pathway to directory containing IQ-Tree output files
  --topobinner          Basic RF-distance based binning of topologies (RF==0) - Provide Tree Viewer input with empty TopologyID column
  --phybin_external     Combine PhyBin output with Tree Viewer file from the IQ-Tree or IQ-Tree_external stages

Trimal arguments (--trimal):
  --trimal_gap_threshold
                        Minimum percent of valid bases per-sequence in alignment required [0 - 1.0] (default: 0.9)
  --trimal_min_seq_len
                        Minimum sequence length for a window to be retained (default: 1kb)
  --trimal_drop_windows
                        Drop windows missing sample sequence (default: True)

Pairwise Estimator arguments (--pw_estimator):
  --pwe_percent_chrom   Percentage of each chromosome to randomly draw windows (default: 0.1)

Pairwise Filter arguments (--pw_filter):
  --pw_subwindow_size   Sub-window size for pairwise filtration [bp/kb/mb] (default: 100bp)
  --pw_step             Step size for pairwise filtration [bp/kb/mb] (default: 10bp)
  --pw_min_seq_len      Minimum sequence length for pairwise filtration [bp/kb/mb] (default: 1000bp)
  --pw_max_pdist        Maximum p-distance per sub-window allowed before sub-window sequence is masked. This only affects the sample with the elevated signal. (default: 0.15)
  --pw_zscore_cutoff    Maximum number of standard deviations away from mean for a windows sequence to be kept in analysis (default: 2)
  --pw_seq_coverage     Minimum pairwaise coverage - window with any sample below cutoff is dropped (default: 0.9)
  --pw_missing_char     Character used to represent masked sequence (default: N)
  --pw_exclude          List of samples to exclude when running pairwise filter (default: Reference)

IQ-Tree arguments (--iqtree):
  --iqtree_model        Model to pass to IQ-Tree (default: GTR*H4)
  --iqtree_bootstrap    Number of bootstrap replicated to pass to IQ-Tree. IQ-Tree requires a minimum of 1000 replicates. (default: 1000)
  --iqtree_cpu_cores    Number of cores to give each IQ-Tree run - recommended to use 3-5 cores max (default: 4)
  --tv_file_name        Name of Tree Viewer input file produced at the end of the Tree Viewer pipeline (default: TreeViewer_input_file.xlsx)

Topobinner arguments (--topobinner):
  --tb_rooted_trees     True if input trees are rooted, False if they are unrooted (default: False)

Tree Viewer pipeline tools:
  --tv_config_template  Output a blank TreeViewer pipeline configuration file
  --parse_treeviewer    Creates individual files for each chromosome of a Tree Viewer input file - (best for large genomes with small window sizes)
  --rootTV [ ...]       Root Newick trees in Tree Viewer input file. List multiple samples with a space between them. (i.e., cat1 cat2 cat3)

p-Distance Tracer tools (--pdistance):
  --pdistance           Calculate p-distance for one or more multiple-sequence alignment fasta files.
  --pdist_filename      File name for p-distance calculator. (default: SignalTracer_input.tsv)
  --pdist_threshold     Maximum frequency of missing data per-window. (default: 0.75)
  --pdist_missing_character
                        Missing data character used to filter samples in --pdistance command (default: N)
  --pdist_ignore_missing
                        Ignore sites with missing data when calculating p-distance (default: True)

System Options:
  --log_level           Log level (NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL)
  --cpu                 Number of CPUs to use (default is system max)
```

## Input Files

All of the pipelines for THExBuilder start with multiple-sequence alignment fasta files, organized within a single directory. Each multiple-sequence alignment fasta file ideally represents a single chromosome, but scaffold level assemblies can also be used with the caveat that there is no way to link scaffolds into pseudo-chromosomes within the pipeline or application. Input fasta files can be compressed using bgzip or uncompressed.
Below is an example of the directory and file structure of the input multiple-sequence fasta files.

```
MultiAlignmentDir/  
| -- chr1.fasta  
| -- chr2.fasta  
| -- chX.fasta  
```

Example file: (chr1.fasta)  

```
> Sample 1  
AGTGCTAGC…GTTC  
> Sample 2  
AGTGTCCTA…GCTT  
> Sample 3  
AGTGTCGTA…GCTT   
```

## Tree Viewer Pipeline

The Tree Viewer pipeline is designed to take multi-alignment fasta files, parse them into non-overlapping windows, run the windows through several filtration steps, infer maximum likelihood phylogenies per-window using IQ-Tree2, and then bin identical topologies using _topobinner_ for visualization in the Tree Viewer dashboard. The examples below show how to run each stage with and without a configuration file, providing the user flexibility to run the pipeline with all arguments conveniently located in a single file. It is recommended to use a configuration file as this helps improve reproducibility and also reduces the length of the commands, but the option is left to the user. Since the configuration file has to be built with a specific format, you can generate a blank configuration file by running “thexb --tv_config_template”.

### Start-to-Finish Run

The Tree Viewer pipeline is designed to be flexible, providing the user the ability to run each step independently or continuously one after another. Although the option to run the pipeline start-to-finish is available, it is recommended to run each step individually to inspect and validate the intermediate results. However, if you choose to run it all at once you will need to create a configuration file and pass it to the '-c' argument.

```
thexb --tv_all -c config.ini  
```

### 1. Fasta Windowing

The fasta windowing step takes a directory containing multi-alignment fasta files and returns a directory of subdirectories, labeled by the original sequence names (i.e., chr1.fasta -> chr1/), containing the non-overlapping windowed fasta files.
Note that this step can return thousands of files, so be mindful when opening these directories in a file browser. Fasta files can be uncompressed or compressed using bgzip.

```
thexb --minifastas -c config.ini  
```

or

```
thexb --minifastas \
-i ./chromosome_alignments \
-o ./output \
--window_size 100kb 
```

### 2. Trimal Gap Trimming

The Trimal stage of the pipeline uses the gap-threshold feature to remove columns of a window alignment where the amount of missing data exceeds the provided threshold.
A second threshold (trimal_min_seq_len) is provided to remove trimmed window alignments shorter than the minimum sequence length defined.

```
thexb --trimal -c config.ini
```

or

```
thexb --trimal \
-i ./output/windowed_fastas \
-o ./output \
--trimal_min_seq_len 1000bp \
--trimal_gap_threshold 0.85
```

### 3. Pairwise Distance + Coverage Filter

The pairwise distance and coverage filter walks through each window of the genome and conducts a sliding sub-window analysis with a step size of n-bp to mask local regions of extreme divergence that may be caused by misalignment or undetected paralogy. It also masks entire windows where a single taxon’s base coverage is below the provided threshold. This step was integrated from Foley et al. 2021 (6) and updated to run on Python 3.

```
thexb --pw_filter -c config.ini
```

or

```
thexb --pw_filter \
-i ./output/trimal_filtered_windows \
-o ./output \
--reference ReferenceName \
--pw_subwindow_size 100bp \
--pw_step 10bp \
--pw_min_seq_len 5000bp \
--pw_max_pdist 1.0 \
--pw_zscore_cutoff 2 \
--pw_missing_char "-" \
--pw_seq_coverage 0.8 
```

### 4. IQ-Tree: Maximum Likelihood Phylogeny Inference

IQ-Tree2 is used to infer maximum likelihood phylogenies, providing the basis for the distribution of phylogenetic signal visualized in the Tree Viewer dashboard. The filtered fasta files generated by the pairwise filtration step are passed to IQ-Tree2 and the resulting Newick files are collected and organized into a Tree Viewer input file. There are two approaches built into the Tree Viewer pipeline for the IQ-Tree analysis. The first is to have THExBuilder call IQ-Tree for each filtered window by using the “--iqtree” command, a model of sequence evolution, and number of bootstrap replications. The second approach is to take the filtered fasta files from the pairwise filtration step and run IQ-Tree externally on a different local machine, server, or cluster. If you choose to run IQ-Tree externally, you will need to organize the resulting “.treefile” Newick trees into a single directory and pass the directory to the “-iqtree_external” command rather than the “--iqtree” command. You can also optionally create sub-directories per chromosome and organize the “.treefile” by chromosome and pass it to the “--iqtree_external” command, but it is not required.

_NOTE: minimum required bootstrap replicates is 1000_

```
thexb --iqtree -c config.ini
```

or

```
thexb --iqtree \
-i ./output/pairwise_filtered_windows/ \
-o ./output/ \
--iqtree_model GTR*H4 \
--iqtree_bootstrap 1000 \
--tv_file_name TVinput.xlsx
```

### 5. Topobinner

Topobinner is a tool used to organize and label identical tree topologies based on RF-distance. Trees are treated as equal topologies when their RF-distance is equal to 0 (9-10). The trees are binned and then labeled based on their whole genome frequency, meaning the most frequent topology in the genome will be labeled Tree1, then Tree 2 for the second most frequent topology, and so on. Once this step is completed, an updated Tree Viewer file with binned topologies will be generated and placed in the defined output directory. It should also be noted that users can also pass a Tree Viewer input file produced by the File Pruning export option (see Tree Viewer export options section for more details) in Tree Viewer to bin new trees that have been pruned into a subset from a larger Tree Viewer input file. No special steps need to be taken to run this, simply provide the Tree Viewer input file as you would normally.  

```
thexb --topobinner -c config.ini
```

or

```
thexb --topobinner \
-i ./output/TVinput.xlsx \
-o ./output/
```

### Topobinner Alternative

An alternative to using topobinner is PhyBin (<https://github.com/rrnewton/PhyBin>). PhyBin is not hosted on conda, so you would need to download the program from GitHub and install it externally.
It is important that you use PhyBin's "--bin" option and not it’s clustering algorithm. When you run PhyBin, it will generate binned results in a new output directory named either phybin_output or whatever you provide as an output location. Passing the output directory generated by PhyBin and a Tree Viewer file with the TopologyID column blank to THExBuilder's "--phybin_external" command will produce a complete Tree Viewer input file just like topobinner.

```
thexb --phybin_external \
-i phybin_output_directory/ \
--tv_file_name TVinput.xlsx
```

## Additional Tree Viewer Pipeline Tools

### Root trees in a Tree Viewer file

Tree Viewer files with unrooted trees can quickly be rooted to an outgroup by using the --rootTV option in THExBuilder. Multiple outgroup taxa can be defined by providing the names of samples in a space delimited list after the --rootTV argument.

```
thexb --rootTV sample1 sample2 \
-i TVinput.xlsx
-o ./output/
```

## Signal Tracer p-Distance Pipeline

The Signal Tracer pipeline is a single script that calculates raw, uncorrected p-distance between a reference and other taxon of one or more multi-alignment fasta files. There are two ways to run the pipeline, you can provide a single multi-alignment fasta file or a directory of multi-alignment fasta files. The directory of files is typically the input directory used for the entire Tree Viewer pipeline that is ideally a multi-alignment fasta file for each chromosome. The script enables users to include or exclude missing data in the p-distance calculation by providing the --pdist_ignore_missing argument, and also the --pdist_missing_character if missing data is represented by something other than "N".

```
thexb --pdistance \
-i ./chromosome_alignments \
-o ./output/ \
-r Tiger \
-w 100kb \
--pdist_threshold 0.75
```

# References

1. Edelman, N. B. et al. Genomic architecture and introgression shape a butterfly radiation. Science (80-. ). 366, 594–599 (2019).

2. Li, G., Figueiró, H. V, Eizirik, E., Murphy, W. J. & Yoder, A. Recombination-Aware Phylogenomics Reveals the Structured Genomic Landscape of Hybridizing Cat Species. Mol. Biol. Evol. 36, 2111–2126 (2019).

3. Nelson, T. C. et al. Ancient and recent introgression shape the evolutionary history of pollinator adaptation and speciation in a model monkeyflower radiation (Mimulus section Erythranthe). PLOS Genet. 17, e1009095 (2021).

4. Small, S. T. et al. Radiation with reticulation marks the origin of a major malaria vector. Proc. Natl. Acad. Sci. 202018142 (2020) doi:10.1073/pnas.2018142117.

5. trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.

6. Foley et al. Zoonomia Phylogenetic Analysis (unpublished)

7. L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274.

8. S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2019) GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol., in press. <https://doi.org/10.1093/sysbio/syz051>

9. Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. Math. Biosci. 53, 131–147 (1981).

10. ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046

11. Bredemeyer, K. R. et al. Ultracontinuous Single Haplotype Genome Assemblies for the Domestic Cat (Felis catus) and Asian Leopard Cat (Prionailurus bengalensis). J. Hered. 112, 165–173 (2021).
