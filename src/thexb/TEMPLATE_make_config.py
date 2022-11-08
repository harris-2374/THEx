def config_template():
    file_contents = """[General]
multi_alignment_dir = input/dir/with/chromosomes/
TreeViewer_file_name = TreeViewer_input_file.xlsx
outdir = /example/output/directory

[Fasta Windower]
# Give as 100bp/kb/mb
window_size = 100bp

[Trimal]
gap_threshold = 0.9
minimum_seq_length = 50bp
drop_windows = True

[Pairwise Estimator]
percent_of_chromosome_to_run = 0.1

[Pairwise Filter]
reference_name = REFERENCE
filter_window_size = 100bp
step = 10bp
min_seq_len = 1000bp
max_pDistance_cutoff = 0.024
Zscore = 2
pairwise_coverage_cutoff = 0.9
exclude_list = REFERENCE
missing_char = N

[IQ-TREE]
model = GTR*H4
bootstrap = 1000
cores_per_job = AUTO

[Topobinner]
rooted_trees = N

[Logging]
level = INFO

[Processing]
multiprocess = 8
    """
    with open("config_template.ini", 'w') as oh:
        oh.write(file_contents)
    return