class HelpDesc():
    def iqtree_external(self):
        msg = "Create base TreeViewer input file with externally run IQ-Tree results. Use (-i) to provide pathway to directory containing IQ-Tree output files"
        return msg

    def iqtree(self):
        msg="Run windowed fasta files through IQTree"
        return msg

    def tv_all(self):
        msg="Run all steps of the Tree Viewer pipeline sequentially"
        return msg

    def fasta_windowing(self):
        msg="Fasta windowing step"
        return msg

    def trimal(self):
        msg="Run windowed fasta files through Trimal"
        return msg

    def pw_estimator(self):
        msg="Additional tool that provides basic information to help choose initial settings for the pairwise filtration step"
        return msg

    def pw_filter(self):
        msg="Sliding window approach to masking segments of windows where a sample's divergence is above the provided z-score threshold as well as dropping windows that contain samples with coverage below given threshold."
        return msg

    def topobinner(self):
        msg="Basic RF-distance based binning of topologies (RF==0) - Provide Tree Viewer input with empty TopologyID column"
        return msg
    
    def phybin(self):
        msg="Combine PhyBin output with Tree Viewer file from the IQ-Tree or IQ-Tree_external stages"
        return msg

    def pdistance(self):
        msg="Calculate p-distance for one or more multiple-sequence alignment fasta files."
        return msg
    
    def pdistance_threshold(self):
        msg = "Maximum frequency of missing data per-window. (default: 0.75)"
        return msg

    def pdistance_filename(self):
        msg = "File name for p-distance calculator. (default: SignalTracer_input.tsv)"
        return msg

    def reference(self):
        msg="Reference sample for the --pdistance, --pw_estimator, and --pw_filter commands."
        return msg

    def parse_treeviewer(self):
        msg="Creates individual files for each chromosome of a Tree Viewer input file - (best for large genomes with small window sizes)"
        return msg

    def config_template(self):
        msg="Output a blank TreeViewer pipeline configuration file"
        return msg

    def update_freq(self):
        msg="Number of files to run before giving update"
        return msg

    def config(self):
        msg="Pathway to configuration file for Tree Viewer pipeline (*recommended approach*)"
        return msg

    def input(self):
        msg="Input file path (Cannot be used when running Tree Viewer pipeline with configuration file)"
        return msg

    def output(self):
        msg="Output directory (Cannot be used when running Tree Viewer pipeline with configuration file)"
        return msg

    def log_level(self):
        msg='Log level (NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL)'
        return msg

    def cpu(self):
        msg='Number of CPUs to use (default is system max)'
        return msg

    def window_size(self):
        msg='Window size [bp/kb/mb/gb] (default: 100kb)'
        return msg

    def trimal_gap_threshold(self):
        msg='Minimum percent of valid bases per-sequence in alignment required [0 - 1.0] (default: 0.9)'
        return msg

    def trimal_minSeqLen(self):
        msg='Minimum sequence length for a window to be retained (default: 1kb)'
        return msg

    def pw_window_size(self):
        msg='Sub-window size for pairwise filtration [bp/kb/mb] (default: 100bp)'
        return msg
    
    def pw_step(self):
        msg='Step size for pairwise filtration [bp/kb/mb] (default: 10bp)'
        return msg

    def pw_min_seq_len(self):
        msg='Minimum sequence length for pairwise filtration [bp/kb/mb] (default: 1000bp)'
        return msg

    def pw_max_pdist(self):
        msg='Maximum p-distance per sub-window allowed before sub-window sequence is masked. This only affects the sample with the elevated signal. (default: 0.15)'
        return msg
    
    def pw_zscore_cutoff(self):
        msg='Maximum number of standard deviations away from mean for a windows sequence to be kept in analysis (default: 2)'
        return msg

    def pw_seq_coverage(self):
        msg='Minimum pairwaise coverage - window with any sample below cutoff is dropped (default: 0.9)'
        return msg

    def pw_missing_char(self):
        msg='Character used to represent masked sequence (default: N)'
        return msg
    
    def pw_exclude(self):
        msg='List of samples to exclude when running pairwise filter (default: Reference)'
        return msg

    def pwe_percent_chrom(self):
        msg='Percentage of each chromosome to randomly draw windows (default: 0.1)'
        return msg

    def iqtree_model(self):
        msg='Model to pass to IQ-Tree (default: GTR*H4)'
        return msg

    def iqtree_bootstrap(self):
        msg='Number of bootstrap replicated to pass to IQ-Tree. IQ-Tree requires a minimum of 1000 replicates. (default: 1000)'
        return msg
    
    def iqtree_cpu_cores(self):
        msg='Number of cores to give each IQ-Tree run - recommended to use 3-5 cores max (default: 4)'
        return msg

    def tb_rooted_trees(self):
        msg='True if input trees are rooted, False if they are unrooted (default: False)'
        return msg
    
    def tb_outgroup(self):
        msg='Sample(s) to set as outgroup (default: Will not root trees)'
        return msg

    def parsimony(self):
        msg='***EXPERIMENTAL*** Analyze number of parsimony informative sites per window across the genome - (-i = THExBuilderOutput/windowed_fastas)'
        return msg

    def range(self):
        msg="Region to extract chr:start:stop - (parsimony informative site tool)"
        return msg

    def missing_character(self):
        msg = "Missing data character used to filter samples in --pdistance command (default: N)"
        return msg

    def rootTV(self):
        msg = "Root Newick trees in Tree Viewer input file. List multiple samples with a space between them. (i.e., cat1 cat2 cat3)"
        return msg

    def tv_file_name(self):
        msg = "Name of Tree Viewer input file produced at the end of the Tree Viewer pipeline (default: TreeViewer_input_file.xlsx)"
        return msg

    def trimal_dropwindows(self):
        msg="Drop windows missing sample sequence (default: True)"
        return msg

    def ignore_missing(self):
        msg="Ignore sites with missing data when calculating p-distance (default: True)"
        return msg