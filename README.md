# Tree House Explorer
 Tree House Explorer (thex) is a novel phylogenomic genome browser that host a growing collection of interactive dashboards. It was developed using [Plotly's](https://plotly.com/) open source framework [Dash](https://dash.plotly.com/) for building extremely powerful and customizable analytical web applications.


## THEx on Windows
Two of the main dependencies required for the THExBuilder's Tree Viewer pipeline are not supported on Windows. Due to this limitation, THEx for Windows is only deployed with the main THEx application.

## Browser Support
THEx is supported by most major broswers (Google Chrome, Safari, Firefox, etc.) except Internet Explorer.  


#  Conda Installation

Steps:
1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
    * Windows 10 installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html)
    * Mac OS installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/macos.html)  
    * Linux installation [documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html)
2. Add required channels
    > $ conda config --add channels bioconda
    >
    > $ conda config --add channels conda-forge
3. Create a new conda environment
    > $ conda create -n env_name
4. Activate conda environment
    > $ conda activate env_name
5. Install THEx
    > $ conda install thex -c ajharris_2374  

   
If you run into issues getting miniconda installed and working, there are numerous videos on YouTube that will walk you through the installation process. 
***
# Tree House Explorer Start-up:
Once you have successfully installed Tree House Explorer, it is super simple to get the application up and running in your browser. All you need to do is type the command shown below... That's it!

    $ thex


Once you hit enter, it may take a moment to start up but in the console you should see "Tree House Explorer running on http://127.0.0.1:8050/" pop up. <br>

You can click or copy/paste this address into your browser and the homepage should load up. 

Now that you are on the homepage for THEx, you are ready to begin exploring! If it is your first time using THEx, it's reccomended to take a moment and look over the Docs to get a more detailed breakdown of each dashboard, as well as THEx's built-in command line suite _THExBuilder_. 


# THExBuilder
## Purpose
THExBuilder is a suite of pipelines and tools designed to generate, modify, and analyze THEx input files. It provides a customized solution for generating Tree Viewer input files. It also contains a single command pipeline for calculating raw p-distance values from multiple-sequence fasta files in non-overlapping sliding windows for Signal Tracer. THExBuilder was designed to be flexible and scalable, taking advantage of multiprocessing and partitioned workflows. For example, the Tree Viewer Pipeline is designed to be run either with a configuration file or by command line arguments, as well as being able to run the pipeline step-by-step or from start-to-finish with a single command. This provides the user the freedom and flexibility to tailor their runs and track the progress through pipelines with ease.
## Usage
In addition to being a novel interactive genome browser, THEx comes with a built-in command line suite that contains tools and pipelines that help simplify the process of building and manipulating input files for the dashboards. There are two main pipelines, one for Tree Viewer and one for Signal Tracer. To use THExBuilder, activate the conda environment that THEx is installed within and enter the command _thexb_ to begin. This command will bring up the help section, providing an overview of all the options and available tools. 

THExBuilder's output is contained in a single THExBuilderOutput directory. By default, the output directory is written to the current working directory, so if you do not provide an output location be sure to use the same working directory as you move through the pipelines. Within the output directory you will find the intermediate output files for each stage of the different pipelines as well as log files to help track what has been performed, to what files, in what order, and with what input parameters. To use THExBuilder, activate the conda environment that THEx is installed in and run _thexb_ to bring up the help menu.

## Input Files
All of the pipelines for THExBuilder start with multiple-sequence alignment fasta files, organized within a single directory. Each fasta file ideally represents the multiple-sequence alignment of a single chromosome, but scaffold level assemblies can be used too with the caveat that each sequence will be treated as independent sequences and there is no way to link scaffolds into pseudo-chromosomes within the pipeline or application. Below is a basic example of the directory structure and file structure of the multiple-sequence fasta file.
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
The Tree Viewer pipeline is designed to take multi-alignment fasta files, run them through several filtration steps, and then through IQ-Tree to generate per-window maximum likelihood phylogenies. These phylogenies along with their chromosome-window coordinates are put together in a Tree Viewer file where the topologies can then be binned and labeled for viewing in Tree Viewer. The examples below show how to run each stage with and without a configuration file, providing the user flexibility to run the pipeline with all arguments conveniently located in a single file. It is recommended to use a configuration file as this helps improve reproducibility and also reduces the length of the commands, but the option is left to the user. Since the configuration file has to be built with a specific format, you can generate a blank configuration file running “thexb --tv_config_template”.
### - Start-to-Finish Run -
The Tree Viewer pipeline is designed to be flexible, providing the user complete control over how the pipeline is ran. Although not recommended, one can run the entire Tree Viewer pipeline start to finish by providing a configuration file and calling the “--tv_all” argument. This will run through each of the pipeline steps described below, but will not stop after each step is completed. Although the option to run the pipeline start-to-finish is available, it is recommended to run each step individually to inspect and validate the intermediate results. 

Command Example:
```
$ thexb –tv_all -c config.ini  
```

### 1. Fasta Windowing
The fasta windowing stage will take a directory containing multi-alignment fastas (ideally whole chromosomes) and will return a directory of subdirectories, labeled by the original sequence names (i.e., chromosome names), where each subdirectory contains the windowed fasta files. Note that this step can return thousands of files when using a small window size (i.e., 5kb for a 2.5Gb genome), so be mindful when opening these directories in a file browser.

Command Example:
```
$ thexb --minifastas -c config.ini  
```
or  
```
$ thexb --minifastas --window_size 10kb -i dir_of_multi_alignment_files/   
```

Descriptions:
```
--window_size = Window size [bp/kb/mb] (default: 100kb)
-i = Input file path (Cannot be used when running Tree Viewer pipeline with configuration file)
```

### 2. Trimal Gap Trimming 
The Trimal stage of the pipeline uses the gap-threshold feature that removes regions of a window where a taxon has missing data that exceeds the provided threshold for allowed gaps. A second threshold is provided to remove resulting sequences that do not meet a minimum sequence length.
Command Example:
```
$ thexb --trimal -c config.ini
```
or
```
$ thexb --trimal –trima_gapl_threshold 0.1 --trimal_min_seq_len 1kb -i THExBuilderOutput/windowed_fastas
```

Descriptions:
```
--trimal_gap_threshold = Percentage of alignment with no gaps (i.e., removes sequence with gaps in 10% or more of the sequences
--trimal_min_seq_len = Minimum sequence length to retain window in study. If above, window is excluded
-i = Input file path (Cannot be used when running Tree Viewer pipeline with configuration file)
```

### 3. Pairwise Distance + Coverage Filter 
The pairwise distance and coverage filter walks through each window of the genome and conducts a sliding sub-window analysis with a step size of n-bp to mask local regions of extreme divergence that may be caused by misalignment or undetected paralogy. It also masks entire windows where a single taxon’s base coverage (valid bases include A, T, G, and C) is below the provided threshold. This step was integrated from Foley et al. 2021 (6) and updated to run on Python 3.
Command Example:
```
$ thexb --pw_filter -c config.ini
```
or
```
$ thexb --pw_filter --pw_subwindow_size 100bp --pw_step 10bp --pw_min_seq_len 1000bp --pw_min_pdist 0.15 --pw_zscore_cutoff 2 --pw_seq_coverage 0.9 -i THExBuilderOutput/trimal_filtered_windows/	
```
Descriptions:
```
--pw_subwindow_size = Sub-window size (default: 100bp)
--pw_step 10bp = Sub-window step size (default: 10bp)
 --pw_min_seq_len = Minimum length of sequence in a given window to be retained (default: 1000bp)
--pw_max_pdist = Maximum p-distance for a given taxon to be considered non-extreme. Any sequence above this value in a given sub-window will be masked. Note this value is dataset and taxon dependent! (default: 0.15)
--pw_zscore_cutoff = Maximum number of standard deviations away from the mean (calculated as z-score) allowed before entire window is masked. (default: 2)
--pw_seq_coverage = Minimum valid, non-missing sequence coverage. Individual sequences below this value are masked. (default: 0.9)
-i = Pathway to directory containing the output from the Trimal stage 
(i.e., THExBuilderOutput/trimal_filtered_windows/)
```
### 4. IQ-Tree: Maximum Likelihood Phylogeny Inference 
IQ-Tree is used to infer maximum likelihood phylogenies providing the basis for the distribution of phylogenetic signal visualized in Tree Viewer. The filtered fasta files generated by the pairwise filtration step is passed to IQ-Tree and the resulting Newick files are collected and organized into the basis of the Tree Viewer input file. IQ-Tree creates several output files per-window, but we are only interested in the trees so the other data can be saved or discarded.
There are two approaches built into the Tree Viewer pipeline for the IQ-Tree analysis. The first is to have THExBuilder call IQ-Tree for each filtered window by using the “--iqtree” command, a model of sequence evolution, and number of bootstrap replications. The second approach is to take the filtered fasta files from the pairwise filtration step and run IQ-Tree externally on a different local machine, server, or cluster. If you choose to run IQ-Tree externally, you will need to organize the resulting “.treefile” Newick trees into a single directory and pass the directory to the “-iqtree_external” command rather than the “--iqtree” command. You can also optionally create sub-directories per chromosome and organize the “.treefile” by chromosome and pass it to the “--iqtree_external” command, but it is not required.
Command Example:
```
$ thexb --iqtree -c config.ini
```
or
```
$ thexb --iqtree --iqtree_model “GTR*H4” --iqtree_bootstrap 1000 -i THExBuilderOutput /pairwise_filtered_windows/
```
Descriptions:
```
--iqtree_model = Any valid model of nucleotide sequence evolution supported by IQ-Tree
--iqtree_bootstrap = Number of bootstrap replicates to run
-i = Pathway to directory containing the output from the Pairwise Distance + Coverage Filter step 
```

### 5. Topobinner
Topobinner is a tool used to organize and label equal tree topologies based on RF-distance. Trees are treated as equal topologies when their RF-distance is equal to 0 (9-10). This is what allows you to visualize discordant topologies across the genome and identify regions of interest. The trees are binned and then labeled by their whole genome frequency, meaning the most frequent topology in the genome will be labeled Tree1, then Tree 2 for the second most frequent topology, and so on. Once this step is completed, an updated Tree Viewer file with binned topologies will be generated and placed in the THExBuilderOutput directory ready for visualization in Tree Viewer. It should also be noted that users can also pass a Tree Viewer input file produced by the File Pruning export option (see Tree Viewer export options section for more details) in Tree Viewer to bin new trees that have been pruned into a subset from a larger Tree Viewer input file. No special steps need to be taken to run this, simply provide the Tree Viewer input file as you would normally.
Command Example:
```
$ thexb --topobinner -c config.ini
```
or
```
$ thexb --topobinner -i TreeViewer_file_with_blank_TopologyID_col.xlsx
```

Descriptions:
```
-i = Tree Viewer input file from IQ-Tree step or File Pruning export option with TopologyID column blank
```
#### TOPOBINNER ALTERNATIVE:
An alternative to using topobinner is PhyBin (https://github.com/rrnewton/PhyBin). PhyBin is not hosted on conda, so you would need to download the program from GitHub and install it externally. 
It is important that you use PhyBin's "--bin" option and not it’s clustering algorithm. When you run PhyBin, it will generate binned results in a new output directory named either phybin_output or whatever you provide as an output location. Passing the output directory generated by PhyBin and a Tree Viewer file with the TopologyID column blank to THExBuilder's "--phybin_external" command will produce a complete Tree Viewer input file just like topobinner.
Command Example:
```
$ thexb --phybin_external -i pathway_to_phybin_output_directory/ --tv_file_name ./THExBuilderOutput/TreeViewer_input_file.xlsx
```
### - Signal Tracer Pipeline -
The Signal Tracer pipeline is a single script that calculates raw, uncorrected p-distance between a reference and other taxon of one or more multi-alignment fasta files. There are two ways to run the pipeline, you can provide a single multi-alignment fasta file or a directory of multi-alignment fasta files. The directory of files is typically the input directory used for the entire Tree Viewer pipeline that is ideally a multi-alignment fasta file for each chromosome. Since some chromosomes can be quite large, this process may take several hours to run, if not a day or two depending on the system configuration being used (i.e., processor, amount of memory, etc.). The only two additional arguments that are needed to run the pipeline are the window size and missing data threshold. The missing data threshold is a value between 0.0 and 1.0 that indicates the maximum amount of missing data (gaps or masked sequence) allowed for a single taxon, and if that threshold is surpassed, then the values for all taxa in the window are set to NULL. The output from this pipeline is formatted to be directly uploaded into the Signal Tracer dashboard for quick visualization!
Command Example:
```
$ thexb --pdistance --pdist_threshold 0.75 --window_size 100kb
```
Descriptions:
```
--pdist_threshold = Maximum frequency of missing data allowed in a single taxon before the entire window is masked. (value: 0.0-1.0)
--window_size = Window size [bp/kb/mb] (default: 100kb)
```

# <br>References:

1.	Edelman, N. B. et al. Genomic architecture and introgression shape a butterfly radiation. Science (80-. ). 366, 594–599 (2019).

2.	Li, G., Figueiró, H. V, Eizirik, E., Murphy, W. J. & Yoder, A. Recombination-Aware Phylogenomics Reveals the Structured Genomic Landscape of Hybridizing Cat Species. Mol. Biol. Evol. 36, 2111–2126 (2019).

3.	Nelson, T. C. et al. Ancient and recent introgression shape the evolutionary history of pollinator adaptation and speciation in a model monkeyflower radiation (Mimulus section Erythranthe). PLOS Genet. 17, e1009095 (2021).

4.	Small, S. T. et al. Radiation with reticulation marks the origin of a major malaria vector. Proc. Natl. Acad. Sci. 202018142 (2020) doi:10.1073/pnas.2018142117.

5.	trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. Salvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon. Bioinformatics 2009 25: 1972-1973.

6.	Foley et al. Zoonomia Phylogenetic Analysis (unpublished)

7.	L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. 

8.	S.M. Crotty, B.Q. Minh, N.G. Bean, B.R. Holland, J. Tuke, L.S. Jermiin, A. von Haeseler (2019) GHOST: Recovering historical signal from heterotachously-evolved sequence alignments. Syst. Biol., in press. https://doi.org/10.1093/sysbio/syz051

9.	Robinson, D. F. & Foulds, L. R. Comparison of phylogenetic trees. Math. Biosci. 53, 131–147 (1981).

10.	ETE 3: Reconstruction, analysis and visualization of phylogenomic data. Jaime Huerta-Cepas, Francois Serra and Peer Bork. Mol Biol Evol 2016; doi: 10.1093/molbev/msw046

11.	Bredemeyer, K. R. et al. Ultracontinuous Single Haplotype Genome Assemblies for the Domestic Cat (Felis catus) and Asian Leopard Cat (Prionailurus bengalensis). J. Hered. 112, 165–173 (2021).

