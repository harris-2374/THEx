# -*- coding: utf-8 -*-
"""
Author: Andrew Harris
"""
# --- Native Python imports ---
import argparse
import configparser
import logging
import multiprocessing
import os
import shutil
import sys
from pathlib import Path
# --- Toolkit util imports ---
from thexb.UTIL_converters import convert_window_size_to_int
from thexb.UTIL_parsers import break_alignment_to_chroms
from thexb.UTIL_help_descriptions import HelpDesc
# --- Toolkit pipeline stage imports ---
from thexb.STAGE_minifastas import fasta_windower
from thexb.STAGE_iqtree import iq_tree
from thexb.STAGE_iqtree_external import iq_tree_external
from thexb.STAGE_pairwise_estimator import pairwise_estimator
from thexb.STAGE_pairwise_filter import pairwise_filter
from thexb.STAGE_pdistance_calculator import pdistance_calculator
# from STAGE_percent_missing import percent_missing
from thexb.STAGE_trimal import trimal
from thexb.STAGE_topobinner import topobinner
from thexb.STAGE_phybin import phybin
# --- Toolkit additional tools ----
from thexb.TOOL_parse_treeviewer_per_chromosome import parse_treeviewer_per_chromosome
from thexb.TOOL_parsimony_informative_sites import parsimony_informative_sites
from thexb.TOOL_root_TreeViewer_file import root_TreeViewer_file
from thexb.TOOL_tree_viewer_input_stats import tv_input_stats
# --- Toolkit Template imports ---
from thexb.TEMPLATE_make_config import config_template


############################### Set up logger #################################
logger = logging.getLogger(__name__)
def set_logger(WORKING_DIR, LOG_LEVEL):
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(WORKING_DIR / 'logs/' / 'THExBuilder.log')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    logger.setLevel(LOG_LEVEL)
    return logger

############################# Custom Exceptions ###############################
class Error(Exception):
    """Base class for other exceptions"""
    pass
class ConfigNotFound(Error):
    """Raise when config file is required but not present"""
    pass
class InputNotFound(Error):
    """Raise when input file is required but not present"""
    pass
class InputNotProvided(Error):
    """Raise when not input is provided"""
    pass
class InvalidInput(Error):
    """Raise when -i is provided and --tv_all is called"""
    pass
class InvalidIQTREERoot(Error):
    """Raise when IQ-TREE rooted input from -c is not Y or N"""
    pass
class InvalidMultiprocess(Error):
    "Raise when Multiprocess input is not an integer"
    pass
############################## Helper Functions ###############################
def _check_log_level(l):
    valid_types = {'NOTSET': 0, 'DEBUG': 10, 'INFO': 20, 'WARNING': 30, 'ERROR': 40, 'CRITICAL': 50}
    if l in valid_types.keys():
        return valid_types[l]
    else:
        return False


def _check_input(INPUT):
    try:
        if not INPUT:
            raise InputNotFound
        elif not os.path.exists(INPUT):
            raise InputNotFound
        else:
            INPUT = Path(INPUT)
            return
    except InputNotFound:
        print("ERROR: FileNotFound: Pathway provided to (-i) could not be found or is not provided. Please check input pathway and try again.")
        exit(1)

############################### Main Function #################################
def main():
    """ This section of code handles all input data and calls methods. """
    """===================================================================="""
    # --- Argparse Setup ---
    parser = argparse.ArgumentParser(prog="THExBuilder", description='Command line suite of pipelines and tools for Tree House Explorer')
    general = parser.add_argument_group('General Inputs')
    tv_pipeline = parser.add_argument_group('Tree Viewer Pipeline Stages')
    tv_pipeline_opts = parser.add_argument_group('TreeViewer Pipeline Arguments (no config file)')
    tv_additional_tools = parser.add_argument_group('TreeViewer Pipeline Additional Tools')
    pdist_pipeline = parser.add_argument_group('p-Distance Tracer Pipeline')
    sys_options = parser.add_argument_group('System Options')
    # General inputs (ignored when config file is provided)
    general.add_argument(
        '-i',
        '--input',
        action="store",
        help=HelpDesc().input(),
        default=None,
        metavar='',
    )
    general.add_argument(
        '-o',
        '--output',
        type=str,
        action="store",
        help=HelpDesc().output(),
        default=None,
        metavar='',
    )
    general.add_argument(
        '-r',
        '--reference',
        type=str,
        action='store',
        help=HelpDesc().reference(),
        default=None,
        metavar='',
    )
    general.add_argument(
        '-w',
        '--window_size',
        type=str,
        action="store",
        help=HelpDesc().window_size(),
        default='100kb',
        metavar='',
    )
    general.add_argument(
        '-c',
        '--config',
        type=str,
        action='store',
        help=HelpDesc().config(),
        metavar='',
    )
    general.add_argument(
        '--range',
        type=str,
        action='store',
        help=HelpDesc().range(),
        default=None,
        metavar='',
    )
    general.add_argument(
        '--missing_character',
        type=str,
        action='store',
        help=HelpDesc().missing_character(),
        default='N',
        metavar='',
    )
    # Tree Viewer Pipeline arguments
    tv_pipeline.add_argument(
        '--tv_all',
        action="store_true",
        help=HelpDesc().tv_all(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--minifastas',
        action="store_true",
        help=HelpDesc().fasta_windowing(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--trimal',
        action="store_true",
        help=HelpDesc().trimal(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--pw_estimator',
        action="store_true",
        help=HelpDesc().pw_estimator(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--pw_filter',
        action="store_true",
        help=HelpDesc().pw_filter(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--iqtree',
        action="store_true",
        help=HelpDesc().iqtree(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--iqtree_external',
        action="store_true",
        help=HelpDesc().iqtree_external(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--topobinner',
        action="store_true",
        help=HelpDesc().topobinner(),
        default=False,
    )
    tv_pipeline.add_argument(
        '--phybin_external',
        action="store_true",
        help=HelpDesc().phybin(),
        default=False,
    )
    # p-distance calculator
    pdist_pipeline.add_argument(
        '--pdistance',
        action="store_true",
        help=HelpDesc().pdistance(),
        default=False,
    )
    pdist_pipeline.add_argument(
        '--pdist_threshold',
        type=float,
        action="store",
        metavar='',
        help=HelpDesc().pdistance_threshold(),
        default=0.75,
    )
    # Toolkit
    tv_additional_tools.add_argument(
        '--tv_config_template',
        action="store_true",
        help=HelpDesc().config_template(),
        default=False,
    )
    tv_additional_tools.add_argument(
        '--parse_treeviewer',
        action="store_true",
        help=HelpDesc().parse_treeviewer(),
        default=False,
    )
    tv_additional_tools.add_argument(
        '--rootTV',
        action="store",
        nargs="*",  # 0 or more values expected => creates a list
        metavar='',
        type=str,
        help=HelpDesc().rootTV(),
        default=False,
    )
    tv_additional_tools.add_argument(
        '--parsimony_summary',
        action="store_true",
        help=HelpDesc().parsimony(),
        default=False,
    )
    # Trimal
    tv_pipeline_opts.add_argument(
        '--trimal_gap_threshold',
        type=str,
        action="store",
        help=HelpDesc().trimal_gap_threshold(),
        default=0.9,
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--trimal_min_seq_len',
        type=str,
        action="store",
        help=HelpDesc().trimal_minSeqLen(),
        default='1kb',
        metavar='',
    )
    # Pairwise estimator
    tv_pipeline_opts.add_argument(
        '--pwe_percent_chrom',
        type=float,
        action="store",
        help=HelpDesc().pwe_percent_chrom(),
        default=0.1,
        metavar='',
    )
    # Pairwise filter
    tv_pipeline_opts.add_argument(
        '--pw_subwindow_size',
        type=str,
        action="store",
        help=HelpDesc().pw_window_size(),
        default='100bp',
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_step',
        type=str,
        action="store",
        help=HelpDesc().pw_step(),
        default='10bp',
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_min_seq_len',
        type=str,
        action="store",
        help=HelpDesc().pw_min_seq_len(),
        default='1000bp',
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_max_pdist',
        type=float,
        action="store",
        help=HelpDesc().pw_max_pdist(),
        default=0.15,
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_zscore_cutoff',
        type=float,
        action="store",
        help=HelpDesc().pw_zscore_cutoff(),
        default=0.15,
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_seq_coverage',
        type=float,
        action="store",
        help=HelpDesc().pw_seq_coverage(),
        default=0.9,
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--pw_missing_char',
        type=str,
        action="store",
        help=HelpDesc().pw_missing_char(),
        default='N',
        metavar='',
    )
    # IQ-Tree
    tv_pipeline_opts.add_argument(
        '--iqtree_model',
        type=str,
        action="store",
        help=HelpDesc().iqtree_model(),
        default='GTR*H4',
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--iqtree_bootstrap',
        type=int,
        action="store",
        help=HelpDesc().iqtree_bootstrap(),
        default=1000,
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--iqtree_cpu_cores',
        type=str,
        action="store",
        help=HelpDesc().iqtree_cpu_cores(),
        default='AUTO',
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--tb_rooted_trees',
        type=str,
        action="store",
        help=HelpDesc().tb_rooted_trees(),
        default="N",
        choices=["Y", "N"],
        metavar='',
    )
    tv_pipeline_opts.add_argument(
        '--tv_file_name',
        type=str,
        action="store",
        help=HelpDesc().tv_file_name(),
        default="./THExBuilderOutput/TreeViewer_input_file.xlsx",
        metavar='',
    )
    # Dev Tools
    sys_options.add_argument(
        '--log_level',
        type=str,
        action="store",
        help=HelpDesc().log_level(),
        default='INFO',
        metavar='',
    )
    sys_options.add_argument(
        '--cpu',
        type=int,
        action="store",
        help=HelpDesc().cpu(),
        default=os.cpu_count(),
        metavar='',
    )
    args = parser.parse_args()
    """===================================================================="""
    # --- Stages ---
    ALL_STEPS = args.tv_all
    MINIFASTAS = args.minifastas
    TRIMAL = args.trimal
    PW_ESTIMATOR = args.pw_estimator
    PW_FILTER = args.pw_filter
    IQTREE = args.iqtree
    IQTREE_EXTERNAL = args.iqtree_external
    TOPOBIN = args.topobinner
    PHYBIN = args.phybin_external

    # --- General inputs ---
    WINDOW_SIZE = str(args.window_size)
    RANGE = args.range
    MISSING_CHAR = str(args.missing_character)

    # --- Tree Viewer inputs ---
    TRIMAL_GAP_THRESH = float(args.trimal_gap_threshold)
    TRIMAL_MIN_SEQ_LEN = str(args.trimal_min_seq_len)
    
    PW_WINDOW_SIZE = str(args.pw_subwindow_size)
    PW_STEP = str(args.pw_step)
    PW_MIN_SEQ_LEN = str(args.pw_min_seq_len)
    PW_PDIST_CUTOFF = float(args.pw_max_pdist)
    PW_ZSCORE = float(args.pw_zscore_cutoff)
    PW_SEQ_COV = float(args.pw_seq_coverage)
    PW_MISSING_CHAR = str(args.pw_missing_char)
    PW_ESTIMATOR_PERCENT_CHROM = float(args.pwe_percent_chrom)

    IQTREE_MODEL = str(args.iqtree_model)
    IQTREE_BOOTSTRAP = int(args.iqtree_bootstrap)
    IQTREE_CPU_NUM = str(args.iqtree_cpu_cores)

    TOPOBIN_ROOTED = args.tb_rooted_trees
    # --- Tree Viewer Toolkit inputs ---
    PARSE_TREEVIEWER_FILE = args.parse_treeviewer
    PARSIMONY = args.parsimony_summary
    TV_OUTGROUP = args.rootTV

    # --- p-distance ---
    P_DISTANCE = args.pdistance
    REFERENCE = args.reference
    PDIST_THRESHOLD = float(args.pdist_threshold)

    # --- Config + Logging ---
    CONFIG_FILE = args.config
    LOG_LEVEL = args.log_level
    MULTIPROCESS = args.cpu

    # --- Additional Tools ---
    CONFIG_TEMPLATE = args.tv_config_template
        
    # --- Non-pipeline Input + Ouput ---
    INPUT = args.input
    OUTPUT = args.output

    # Set the working directory
    if not OUTPUT:
        WORKING_DIR = Path.cwd() / 'THExBuilderOutput'
    else:
        WORKING_DIR = Path(OUTPUT) / 'THExBuilderOutput'
    WORKING_DIR.mkdir(parents=True, exist_ok=True)
    # --- Set up logging ---
    log_dir = WORKING_DIR / "logs/"
    log_dir.mkdir(parents=True, exist_ok=True)
    # --- Create Template Files ---
    if CONFIG_TEMPLATE:
        config_template()
        exit()
    # If no argument provided, print help message
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    """===================================================================="""
    if CONFIG_FILE:
        # Read config parser + setup input variables
        config = configparser.ConfigParser()
        try:
            if not CONFIG_FILE:
                raise ConfigNotFound
            else:
                config.read(Path(CONFIG_FILE))
        except ConfigNotFound:
            print("ERROR: FileNotFound: Configuration file called using '-c', but file could not be found. Please check input and rerun")
            exit(1)
        
        # General Input Variables
        INPUT = Path(config['Input Files']['multi_alignment_dir'])
        TREEVIEWER_FN = Path(config['Input Files']['TreeViewer_file_name'])
        MULTIPROCESS = config['Processing']['multiprocess']
        try:
            if not MULTIPROCESS.isnumeric():
                raise InvalidMultiprocess
            else:
                MULTIPROCESS = int(MULTIPROCESS)
        except InvalidMultiprocess:
            print("Invalid input for Multiprocess option. Value must an integer 1-n, where n is the total number of cores available")

        WINDOW_SIZE_STR = config['Fasta Windower']['window_size']
        WINDOW_SIZE_INT = convert_window_size_to_int(WINDOW_SIZE_STR)

        # Trimal
        TRIMAL_THRESH = float(config['Trimal']['gap_threshold'])
        TRIMAL_MIN_LENGTH = convert_window_size_to_int(config['Trimal']['minimum_seq_length'])

        # Pairwise Filtering Input Variables
        PW_WINDOW_SIZE = convert_window_size_to_int(config['Pairwise Filter']['filter_window_size'])
        PW_STEP = convert_window_size_to_int(config['Pairwise Filter']['step'])
        PW_REF = config['Pairwise Filter']['reference_name']
        PW_MIN_SEQ_LEN = convert_window_size_to_int(config['Pairwise Filter']['min_seq_len'])
        PW_PDIST_CUTOFF = float(config['Pairwise Filter']['max_pDistance_cutoff'])
        PW_ZSCORE = float(config['Pairwise Filter']['Zscore'])
        PW_PC_CUTOFF = float(config['Pairwise Filter']['pairwise_coverage_cutoff'])
        PW_EXCLUDE_LIST = config['Pairwise Filter']['exclude_list']
        PW_MISSING_CHAR = config['Pairwise Filter']['missing_char']

        # Pairwise Estimator Input Variables
        PW_EST_PERCENT_CHROM = float(config['Pairwise Estimator']['percent_of_chromosome_to_run'])

        # IQ-TREE Input Variables
        IQT_MODEL = str(config['IQ-TREE']['model'])
        IQT_BOOTSTRAP = int(config['IQ-TREE']['bootstrap'])
        IQT_CORES = str(config['IQ-TREE']['cores_per_job'])

        # TopoBin Input Variables - [Y/N]
        TOPOBIN_ROOTED = str(config['Topobinner']['rooted_trees'])
        try:
            if TOPOBIN_ROOTED == "Y":
                pass
            elif TOPOBIN_ROOTED == "N":
                pass
            else:
                raise InvalidIQTREERoot
        except InvalidIQTREERoot:
            print(f"ERROR: Topobinner 'rooted_trees' input can only be Y or N")
            exit(1)

        # TreeViewer output file
        TREEVIEWER_FN = WORKING_DIR / 'TreeViewer_input.xlsx'

        """===================================================================="""
        # Check log level input + return log level int
        try:
            LOG_LEVEL = _check_log_level(LOG_LEVEL)
            assert type(LOG_LEVEL) == int
        except AssertionError:
            raise AssertionError("Invalid log level - Provide NOTSET, DEBUG, INFO, WARNING, ERROR, or CRITICAL")
        logger = set_logger(WORKING_DIR, LOG_LEVEL)
        """===================================================================="""
        # --- Check inputs for errors ---
        try:
            assert INPUT.is_dir()
        except AssertionError:
            logger.warning("Multi-alignment input is a file, please parse genome into chromosomes and provide path to directory containing per-chromosome multi-alignmnet files")
            raise AssertionError("Multi-alignment input is a file, please parse genome into chromosomes and provide path to directory containing per-chromosome multi-alignmnet files")
    elif INPUT:
        # General Input Variables
        try:
            if TOPOBIN:
                pass
            elif ALL_STEPS:
                raise InvalidInput
            elif not INPUT:
                raise InputNotFound
            elif not os.path.exists(INPUT):
                raise InputNotFound
            else:
                INPUT = Path(INPUT)
                pass
        except InputNotFound:
            print("ERROR: FileNotFound: Pathway provided to (-i) could not be found or is not provided. Please check input pathway and try again.")
            exit(1)
        except InvalidInput:
            print("ERROR: InvalidInput: (-i) argument is not a valid input when running --tv_all, use a configuration file (-c).")
            exit(1)
        
        WINDOW_SIZE_STR = WINDOW_SIZE
        WINDOW_SIZE_INT = convert_window_size_to_int(WINDOW_SIZE_STR)
        # Trimal
        TRIMAL_THRESH = TRIMAL_GAP_THRESH
        TRIMAL_MIN_LENGTH = convert_window_size_to_int(TRIMAL_MIN_SEQ_LEN)
        # Pairwise Filtering Input Variables
        PW_WINDOW_SIZE = PW_WINDOW_SIZE
        PW_STEP = PW_STEP
        PW_REF = REFERENCE
        PW_MIN_SEQ_LEN = PW_MIN_SEQ_LEN
        PW_PC_CUTOFF = PW_SEQ_COV
        PW_ZSCORE = PW_ZSCORE
        PW_PDIST_CUTOFF = PW_PDIST_CUTOFF
        PW_EXCLUDE_LIST = REFERENCE
        PW_MISSING_CHAR = PW_MISSING_CHAR
        # Pairwise Estimator Input Variables
        PW_EST_PERCENT_CHROM = PW_ESTIMATOR_PERCENT_CHROM
        # IQ-TREE Input Variables
        IQT_MODEL = IQTREE_MODEL
        IQT_BOOTSTRAP = IQTREE_BOOTSTRAP
        IQT_CORES = IQTREE_CPU_NUM
        # TopoBin Input Variable set above
        # TreeViewer output file
        TREEVIEWER_FN = Path(args.tv_file_name)

        logger = set_logger(WORKING_DIR, LOG_LEVEL)
    elif TOPOBIN:
        TREEVIEWER_FN = Path(args.tv_file_name)
        logger = set_logger(WORKING_DIR, LOG_LEVEL)
    else:
        try:
            if (not CONFIG_FILE) and (not INPUT):
                raise InputNotProvided
        except InputNotProvided:
            print("ERROR: FileNotFound: No input file provided, please provide an input (-i) or configuration file (-c)")
            exit(1)

    try:
        """===================================================================="""
        # --- p-Distance Tracer Pipeline ---
        if P_DISTANCE:
            # Input/output
            pdistance_output_dir = WORKING_DIR / 'p-distance/'
            pdistance_output_dir.mkdir(parents=True, exist_ok=True)
            # Set up logging dir
            log_dir = WORKING_DIR / "logs/"
            log_dir.mkdir(parents=True, exist_ok=True)
            # Check log level input + return log level int
            try:
                LOG_LEVEL = _check_log_level(LOG_LEVEL)
                assert type(LOG_LEVEL) == int
            except AssertionError:
                raise AssertionError("Invalid log level - Provide NOTSET, DEBUG, INFO, WARNING, ERROR, or CRITICAL")

            # Make sure all required input variables are valid
            try:
                if not REFERENCE:
                    raise AssertionError
            except AssertionError:
                print("ERROR: No reference (-r, --reference) sample provided")
                exit()
            try:
                if not REFERENCE:
                    raise AssertionError
            except AssertionError:
                print("ERROR: No reference (-r, --reference) sample provided")
                exit()
            # Log init + run
            logger.info("=======================================")
            logger.info("======== p-Distance Calculator ======== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info("------------------------------------")
            logger.info(f"Input directory: {INPUT.as_posix()}")
            logger.info(f"Output directory: {pdistance_output_dir.as_posix()}")
            logger.info(f"Reference sample: {REFERENCE}")
            logger.info(f"Missing data threshold: {PDIST_THRESHOLD}")
            logger.info(f"Missing data character: {MISSING_CHAR}")
            logger.info(f"Window size: {WINDOW_SIZE_STR}")
            logger.info("------------------------------------")
            pdistance_calculator(INPUT, pdistance_output_dir, PDIST_THRESHOLD, REFERENCE, MISSING_CHAR, WORKING_DIR, WINDOW_SIZE_INT, MULTIPROCESS, LOG_LEVEL)
            pass
        """===================================================================="""
        # Addition TreeViewer Tools
        if PARSIMONY:
            # Check/set logging
            try:
                LOG_LEVEL = _check_log_level(LOG_LEVEL)
                assert type(LOG_LEVEL) == int
            except AssertionError:
                raise AssertionError("Invalid log level - Provide NOTSET, DEBUG, INFO, WARNING, ERROR, or CRITICAL")
    
            _check_input(INPUT)
            logger.info("============================================= ")
            logger.info("======== Parsimony Informative Sites ======== ")
            logger.info("============================================= ")
            parsimony_output_dir = WORKING_DIR / 'parsimony_informative_site_summary/'
            parsimony_output_dir.mkdir(parents=True, exist_ok=True)
            parsimony_informative_sites(parsimony_output_dir, INPUT, RANGE, WORKING_DIR, LOG_LEVEL)
            pass
        
        if TV_OUTGROUP:
            TV_OUTGROUP = list(TV_OUTGROUP)
            _check_input(INPUT)
            logger.info("============================================= ")
            logger.info("=========== Root Tree Viewer File =========== ")
            logger.info("============================================= ")
            updated_TV_dir = WORKING_DIR / 'updated_TreeViewer_files/'
            updated_TV_dir.mkdir(parents=True, exist_ok=True)
            root_TreeViewer_file(INPUT, updated_TV_dir, TV_OUTGROUP, WORKING_DIR, LOG_LEVEL)

        # --- Tree Viewer Pipeline ---
        if ALL_STEPS:
            logger.info("")
            logger.info("=======================================")
            logger.info("============ Full Pipeline ============ ")
            logger.info("=======================================")
            try:
                if not CONFIG_FILE:
                    raise ConfigNotFound

            except ConfigNotFound:
                print("ERROR: FileNotFound: A configuration file was not provided. A configuration file must be provided when running the full pipeline.")
            MINIFASTAS = True
            TRIMAL = True
            PW_FILTER = True
            IQTREE = True
            TOPOBIN = True
            pass
        
        if MINIFASTAS:
            # Input/output
            outdir = WORKING_DIR / "windowed_fastas"
            # Log input parameters
            logger.info("")
            logger.info("=======================================")
            logger.info("======= Fasta Windowing Parser ======== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Input directory: {INPUT.as_posix()}")
            logger.info(f"Output directory: {outdir.as_posix()}")
            logger.info(f"Window size: {WINDOW_SIZE_STR}")
            logger.info("------------------------------------")
            fasta_windower(INPUT, WORKING_DIR, WINDOW_SIZE_STR, WINDOW_SIZE_INT, MULTIPROCESS, LOG_LEVEL)
            pass

        if TRIMAL:
            # Input/output
            windowed_fasta_dir = WORKING_DIR / 'windowed_fastas'
            filtered_outdir = WORKING_DIR / 'trimal_filtered_windows'
            # Log + run
            logger.info("")
            logger.info("=======================================")
            logger.info("=============== Trimal ================ ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Input directory: {windowed_fasta_dir.as_posix()}")
            logger.info(f"Gap threshold: {TRIMAL_THRESH}")
            logger.info(f"Minimum post-trimal sequence length: {TRIMAL_MIN_LENGTH}")
            logger.info("------------------------------------")
            trimal(windowed_fasta_dir, filtered_outdir, WORKING_DIR, TRIMAL_THRESH, TRIMAL_MIN_LENGTH, MULTIPROCESS, LOG_LEVEL)
            pass

        if PW_ESTIMATOR:
            # Input/output
            filtered_indir = WORKING_DIR / 'trimal_filtered_windows'
            # Log input parameters
            logger.info("")
            logger.info("=======================================")
            logger.info("========== Pairwise Estimator ========= ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Input directory: {filtered_indir.as_posix()}")
            logger.info(f"Reference sample: {PW_REF}")
            logger.info(f"Percentage of chromosome sampled: {PW_EST_PERCENT_CHROM}")
            logger.info("------------------------------------")
            pairwise_estimator(filtered_indir, PW_REF, PW_EST_PERCENT_CHROM, WORKING_DIR, LOG_LEVEL)
            pass
        
        if PW_FILTER:
            # This is the input when percent missing step is active
            # filtered_indir = WORKING_DIR / 'percent_missing_filtered_windowed_chroms'
            
            # Set exluded info into list
            PW_EXCLUDE_LIST = list(PW_EXCLUDE_LIST.replace(' ', '').split(','))
            # Input/output      
            filtered_indir = WORKING_DIR / 'trimal_filtered_windows'
            filtered_outdir = WORKING_DIR / 'pairwise_filtered_windows'
            # Log input parameters
            logger.info("")
            logger.info("=======================================")
            logger.info("=========== Pairwise Filter =========== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Input Directory: {filtered_indir}")
            logger.info(f"Output Directory: {filtered_outdir}")
            logger.info(f"Window Size: {PW_WINDOW_SIZE}")
            logger.info(f"Step Size: {PW_STEP}")
            logger.info(f"Pairwise deletion distance frequency cutoff value: {PW_PDIST_CUTOFF}")
            logger.info(f"Reference Sample: {PW_REF}")
            logger.info(f"Minimum Sequence Length: {PW_MIN_SEQ_LEN}")
            logger.info(f"Pairwise Coverage Cutoff: {PW_PC_CUTOFF}")
            logger.info(f"Z-score: {PW_ZSCORE}")
            logger.info(f"Sample Exclusion List: {PW_EXCLUDE_LIST}")
            logger.info("------------------------------------")
            # Call pairwise_filter
            pairwise_filter(
                filtered_indir,
                filtered_outdir,
                WORKING_DIR,
                PW_WINDOW_SIZE,
                PW_STEP,
                PW_PDIST_CUTOFF,
                PW_REF,
                PW_MIN_SEQ_LEN,
                PW_PC_CUTOFF,
                PW_ZSCORE,
                PW_EXCLUDE_LIST,
                PW_MISSING_CHAR,
                MULTIPROCESS,
                LOG_LEVEL,
            )
            pass

        if IQTREE:
            # Input/output
            filtered_indir = WORKING_DIR / 'pairwise_filtered_windows'
            filtered_metadata_outdir = WORKING_DIR / 'IQ-TREE_metadata'
            filtered_tree_outdir = WORKING_DIR / 'IQ-TREE_trees'
            # Erase old output directories if present and create new directories
            try:
                filtered_metadata_outdir.mkdir(parents=True)
                filtered_tree_outdir.mkdir(parents=True)
            except FileExistsError:
                shutil.rmtree(filtered_metadata_outdir)
                shutil.rmtree(filtered_tree_outdir)
                filtered_metadata_outdir.mkdir(parents=True)
                filtered_tree_outdir.mkdir(parents=True)
            logger.info("")
            logger.info("=======================================")
            logger.info("=============== IQ-TREE =============== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Input directory: {filtered_indir.as_posix()}")
            logger.info(f"IQ-Tree metadata output directory: {filtered_metadata_outdir.as_posix()}")
            logger.info(f"IQ-Tree Newick tree output directory: {filtered_tree_outdir.as_posix()}")
            logger.info(f"Tree Viewer file: {TREEVIEWER_FN}")
            logger.info(f"Model: {IQT_MODEL}")
            logger.info(f"Number of bootstraps: {IQT_BOOTSTRAP}")
            logger.info(f"Number of cores per run: {IQT_CORES}")
            logger.info("------------------------------------")
            iq_tree(
                filtered_indir,
                filtered_metadata_outdir,
                filtered_tree_outdir,
                TREEVIEWER_FN,
                WORKING_DIR,
                IQT_MODEL,
                IQT_BOOTSTRAP,
                IQT_CORES,
                MULTIPROCESS,
                LOG_LEVEL,
            )
            pass

        if IQTREE_EXTERNAL:
            try:
                if not INPUT:
                    raise InputNotFound
                elif not os.path.exists(INPUT):
                    raise InputNotFound
                else:
                    INPUT = Path(INPUT)
                    pass
            except InputNotFound:
                print("ERROR: FileNotFound: Pathway provided to (-i) could not be found or is not provided. Please check input pathway and try again.")
                exit(1)
            logger.info("")
            logger.info("=======================================")
            logger.info("=========== IQ-TREE External ========== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"External input directory: {INPUT.as_posix()}")
            logger.info(f"Tree Viewer file: {TREEVIEWER_FN}")
            logger.info("------------------------------------")
            iq_tree_external(
                INPUT,
                TREEVIEWER_FN,
                WORKING_DIR,
                LOG_LEVEL,
            )
            pass

        if TOPOBIN: 
            # Input/Output
            UPDATED_TV_FILENAME = TREEVIEWER_FN.parents[0] / f"{TREEVIEWER_FN.stem}_topobinner_output.xlsx"
            logger.info("")
            logger.info("=======================================")
            logger.info("========== Topology Binning =========== ")
            logger.info("=======================================")
            logger.info("----------Input Parameters----------")
            logger.info(f"Tree Viewer file: {TREEVIEWER_FN}")
            logger.info(f"Topobinned output file: {UPDATED_TV_FILENAME}")
            logger.info(f"Trees rooted?: {TOPOBIN_ROOTED}")
            logger.info("------------------------------------")
            topobinner(
                TREEVIEWER_FN,
                UPDATED_TV_FILENAME,
                TOPOBIN_ROOTED,
                WORKING_DIR,
                MULTIPROCESS,
                LOG_LEVEL,
            )
        
        if PHYBIN:
            try:
                if not INPUT:
                    raise InputNotFound
                elif not os.path.exists(INPUT):
                    raise InputNotFound
                else:
                    INPUT = Path(INPUT)
                    pass
            except InputNotFound:
                print("ERROR: FileNotFound: Pathway provided to (-i) could not be found or is not provided. Please check input pathway and try again.")
                exit(1)
            logger.info("")
            logger.info("=======================================")
            logger.info("============ PhyBin External ============= ")
            logger.info("=======================================")
            phybin(
                INPUT,
                TREEVIEWER_FN,
                WORKING_DIR,
                LOG_LEVEL,
            )

        # --- THExb Non-pipeline tools ---
        if PARSE_TREEVIEWER_FILE:
            if not OUTPUT:
                WORKING_DIR = Path('TreeViewerFilePerChrom/')
                pass
            else:
                WORKING_DIR = Path(OUTPUT) / 'TreeViewerFilePerChrom/'
                pass
            if not INPUT:
                raise FileNotFoundError("Could not find input file or does not exist, check input and rerun.")
            else:
                INPUT = Path(INPUT)
                pass
            # Create output directory
            WORKING_DIR.mkdir(parents=True, exist_ok=True)
            # Run topobinner
            parse_treeviewer_per_chromosome(
                INPUT,
                WORKING_DIR,
            )
    except KeyboardInterrupt:
        logger.info(f"=======================================")
        logger.info(f"======== Job Cancelled by User ========")
        logger.info(f"=======================================")
    except MemoryError:
        logger.info(f"=======================================")
        logger.info(f"====== Memory Error Encountered ======")
        logger.info(f"=======================================")
    return

