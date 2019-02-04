"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import sys
import inspect
import cemba_data
import logging

log = logging.getLogger()

DESCRIPTION = """
yap (yet another pipeline) is a toolkit for single cell methylation sequencing analysis

This toolkit contain functions for 3-stage analysis:
    Stage 1 Preprocessing: Mapping FASTQ, generate single cell ALLC
    Stage 2 Cell Level Analysis: Prepare MCDS dataset for cell based analysis
    Stage 3 Cluster Level Analysis: merge ALLC and ALLC related functions
    Other functions: qsub submitter for SGE QSUB; simulation functions

STAGE 1
    mapping - Actual mapping function
    default-mapping-config - Print the default mapping config
    mapping-qsub - Qsub wrapper for mapping that run with "yap qsub" 
    mapping-summary - Summary mapping output directory after mapping

STAGE 2
    map-to-region - Map ALLC file into a region BED file
    assemble-dataset - Assemble all cell-region BED file into MCDS
    generate-dataset - Wrapper for map-to-region and assemble-dataset that run with "yap qsub"

STAGE 3
    merge-allc - Merge single cell ALLC files into cluster ALLC.
    allc-profile - Generate summary statistics for a ALLC file.
    allc-to-bigwig - ALLC to BIGWIG.
    allc-extract - Extract ALLC file information. 
    
Qsub
    qsub - Qsubmitter for SGE QSUB.

SIMULATION (exp)
    simulate-allc - Simulate single cell ALLC based on given high coverage ALLC.
    simulate-long-reads-coverage - Simulate genome coverage BED file for long reads.
"""
# TODO add DESCRIPTION for structured functional groups
# mapping
# allc
# dataset etc

EPILOG = """
Author: Hanqing Liu, hanliu@salk.edu
"""


class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)


def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)


def qsub_register_subparser(subparser):
    parser = subparser.add_parser('qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="General qsub helper, need to prepare a command dict file")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the work project"
    )

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the work project"
    )

    parser_req.add_argument(
        "--command_file_path",
        type=str,
        required=True,
        nargs='+',
        help="One or space-separated paths of the command.json file."
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        default=30,
        help="Total CPU in qsub list"
    )

    parser_opt.add_argument(
        "--total_mem",
        type=int,
        required=False,
        default=500,
        help="Total MEM in qsub list"
    )
    return


def pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--fastq_dataframe",
        type=str,
        required=True,
        help="Path of fastq dataframe, can be generate with yap fastq_dataframe"
    )

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser_req.add_argument(
        "--config_path",
        type=str,
        required=True,
        default=None,
        help="Pipeline configuration (.ini) file path. "
             "You can use 'yap default-mapping-config' to print out default config can modify it."
    )
    return


def print_default_config_register_subparser(subparser):
    parser = subparser.add_parser('default-mapping-config',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default config of mapping pipeline")

    parser_opt = parser.add_argument_group("Optional inputs")

    parser_opt.add_argument(
        "--out_path",
        type=str,
        required=False,
        help="Path to save the config file"
    )
    return


def batch_pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping-qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Batch mapping pipeline from multiplexed FASTQ file to ALLC file. "
                                       "Return a command.json file for yap qsub to submit on qsub.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--fastq_dataframe",
        type=str,
        required=True,
        help="Path of fastq dataframe, can be generate with yap fastq_dataframe"
    )

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser_req.add_argument(
        "--config_path",
        type=str,
        required=True,
        default=None,
        help="Pipeline configuration (.ini) file path. "
             "You can use 'yap default-mapping-config' to print out default config can modify it."
    )
    return


def generate_dataset_register_subparser(subparser):
    parser = subparser.add_parser('generate-dataset',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Batch map-to-region for a whole list of allc_files and then generate a MCDS. "
                                       "Generate a command.json file for yap qsub.")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--allc_files",
        type=str,
        required=True,
        help="A path pattern contain wildcard or a path to a file where each line is an ALLC file path"
    )

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser_req.add_argument(
        "--region_bed_path",
        type=str,
        required=True,
        default=None,
        nargs='+',
        help="Space separated path list for reference region bed files"
    )

    parser_req.add_argument(
        "--region_name",
        type=str,
        required=True,
        default=None,
        nargs='+',
        help="Space separated name list for reference region bed files, corresponding to region_bed_path"
    )

    parser_req.add_argument(
        "--context_pattern",
        type=str,
        required=True,
        default=None,
        nargs='+',
        help="Space separated mC context pattern list"
    )

    parser_req.add_argument(
        "--genome_size_path",
        type=str,
        required=True,
        default=None,
        help="File path for the chrom.sizes file of reference genome fasta"
    )

    parser_req.add_argument(
        "--max_cov_cutoff",
        type=int,
        required=True,
        default=None,
        help="Maximum base cov, for normal single cell data, recommended value is 2."
    )

    parser_req.add_argument(
        "--dataset_name",
        type=str,
        required=True,
        default=None,
        help="Name of the dataset output"
    )

    return


def assemble_dataset_register_subparser(subparser):
    parser = subparser.add_parser('assemble-dataset',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Assemble dataset after got all region count bed files for a whole dataset.")
    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    # out_dir, dataset_name, cpu
    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="out_dir path that contains all region count bed files, same as map-to-region step"
    )

    parser_req.add_argument(
        "--region_bed_path",
        type=str,
        required=True,
        default=None,
        nargs='+',
        help="Space separated path list for reference region bed files"
    )

    parser_req.add_argument(
        "--region_name",
        type=str,
        required=True,
        default=None,
        nargs='+',
        help="Space separated name list for reference region bed files, corresponding to region_bed_path"
    )

    parser_req.add_argument(
        "--dataset_name",
        type=str,
        required=True,
        help="Name of the dataset"
    )

    parser_opt.add_argument(
        "--thread",
        type=int,
        required=False,
        default=5,
        help="Number of threads to use, 5 is recommend."
    )


def merge_allc_register_subparser(subparser):
    parser = subparser.add_parser('merge-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Merge ALLC files, will be faster if files all have tabix, "
                                       "otherwise use the old .idx files")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_paths",
        type=str,
        required=True,
        nargs='+',
        help="Space separated ALLC paths OR a file contains all ALLC paths rows. "
             "If provide only one path, each row in the file should be a path"
    )

    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Output path for the merged ALLC file"
    )

    parser_req.add_argument(
        "--chrom_size_file",
        type=str,
        required=True,
        help="Output path for the merged ALLC file"
    )

    parser_opt.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=10,
        help="Number of CPUs for merge ALLC, parallel on genome bins level."
    )

    parser_opt.add_argument(
        "--bin_length",
        type=int,
        required=False,
        default=10000000,
        help="Length of each genome bin to parallel, use default is usually fine."
    )


def map_to_region_register_subparser(subparser):
    parser = subparser.add_parser('map-to-region',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Map ALLC file to region BED, "
                                       "get base mC and coverage count for each region.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="One ALLC file path"
    )

    parser_req.add_argument(
        "--out_path_prefix",
        type=str,
        required=True,
        help="Output file prefix"
    )

    parser_req.add_argument(
        "--region_bed_path",
        type=str,
        required=True,
        nargs='+',
        help="Space separated region BED file paths"
    )

    parser_req.add_argument(
        "--region_name",
        type=str,
        required=True,
        nargs='+',
        help="Space separated region set names corresponding to --region_bed_path"
    )

    parser_req.add_argument(
        "--genome_size_path",
        type=str,
        required=True,
        help="UCSC genome size file"
    )

    parser_req.add_argument(
        "--context_pattern",
        type=str,
        required=True,
        nargs='+',
        help="Space separated methylation context pattern, N for ATCG, H for ATC"
    )

    parser_opt.add_argument(
        "--max_cov_cutoff",
        type=int,
        required=False,
        default=None,
        help="Maximum cutoff for coverage in each base, "
             "e.g. 2 for single cell data, None for bulk seq."
    )

    parser_opt.add_argument(
        "--remove_tmp",
        type=bool,
        required=False,
        default=True,
        help="Remove tmp file or not"
    )
    return


def allc_to_bigwig_register_subparser(subparser):
    parser = subparser.add_parser('allc-to-bigwig',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Just a wrapper of methylpy allc-to-bigwig")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="ALLC path"
    )
    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Out path"
    )
    parser_req.add_argument(
        "--chrom_size",
        type=str,
        required=True,
        help="UCSC chrom.sizes format indicating genome size. ALLC Chr not in this file will be removed."
    )

    parser_opt.add_argument(
        "--mc_type",
        type=str,
        required=False,
        default='CGN',
        help="mC context pattern to use"
    )

    return


def simulate_read_genome_cov_register_subparser(subparser):
    parser = subparser.add_parser('simulate-long-reads-coverage',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Given reads number and length distribution, "
                                       "Simulate a sequencing library in bedgraph format. "
                                       "The last column is genome coverage.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--read_number_mean",
        type=float,
        required=True,
        help="Mean for read number (normal) distribution"
    )
    parser_req.add_argument(
        "--read_number_sd",
        type=float,
        required=True,
        help="SD for read number (normal) distribution"
    )

    parser_req.add_argument(
        "--read_length_mean",
        type=float,
        required=True,
        help="Mean for read length (normal) distribution"
    )

    parser_req.add_argument(
        "--read_length_sd",
        type=float,
        required=True,
        help="SD for read length (normal) distribution"
    )

    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Path of the output file"
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path to the chromosome size file. UCSC chrom.sizes format."
    )

    parser_opt.add_argument(
        "--genome_cov",
        type=int,
        required=False,
        default=2,
        help="Number of genome coverage"
    )

    parser_opt.add_argument(
        "--remove_chr",
        type=bool,
        required=False,
        default=False,
        help="whether remove 'chr' character from the chrom column."
    )


def simulate_allc_register_subparser(subparser):
    parser = subparser.add_parser('simulate-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Given read coverage bedgraph and target ALLC file, "
                                       "simulate a pseudo-ALLC file where (mc / cov) of each base follow "
                                       "beta-binomial distribution.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--genome_cov_path",
        type=str,
        required=True,
        help="Path of the genome coverage bedgraph, the last column should be genome coverage (int) for each region. "
             "Similar to those generated by bedtools genomecov with -bg option."
    )

    parser_req.add_argument(
        "--target_allc_path",
        type=str,
        required=True,
        help="Target ALLC path to simulate from"
    )

    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Path of the output ALLC file"
    )

    parser_opt.add_argument(
        "--allc_profile_path",
        type=str,
        required=False,
        default=None,
        help="Path of the ALLC profile, if None, will search {target_allc_path}.profile, "
             "if not found, an error will occur."
    )


def allc_profile_register_subparser(subparser):
    parser = subparser.add_parser('allc-profile',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Given reads number and length distribution, "
                                       "Simulate a sequencing library in bedgraph format. "
                                       "The last column is genome coverage.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="Path of the ALLC file."
    )

    parser_opt.add_argument(
        "--drop_n",
        type=bool,
        required=False,
        default=True,
        help="Whether to drop context that contain N, such as CCN. "
             "This is usually very rare and need to be dropped."
    )

    parser_opt.add_argument(
        "--n_rows",
        type=int,
        required=False,
        default=100000000,
        help="Number of rows to calculate the profile from. "
             "1e8 rows is a sufficient number to get pretty precise assumption."
    )

    parser_opt.add_argument(
        "--out_path",
        type=str,
        required=False,
        default=None,
        help="Path of the output file. If None, use {allc_path}.profile"
    )


def extract_context_allc_register_subparser(subparser):
    parser = subparser.add_parser('allc-extract',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Extract sub-ALLC for certain mc context, such as mCG")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="Path of the ALLC file."
    )

    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Path of the output ALLC file."
    )

    parser_opt.add_argument(
        "--merge_strand",
        type=bool,
        required=False,
        default=True,
        help="Whether to merge the continuous +/- strand mC and cov, only valid for mCG"
    )

    parser_opt.add_argument(
        "--mc_context",
        type=str,
        required=False,
        default=['CGN'],
        nargs='+',
        help="space separated mC contexts to extract"
    )


def standardize_allc_register_subparser(subparser):
    parser = subparser.add_parser('allc-standardize',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Standardize ALLC files in a whole directory by: "
                                       "1. use bgzip; "
                                       "2. tabix ALLC file; "
                                       "3. remove header line because ALLC columns are always fixed; "
                                       "4. fix 'chr' problem in first column. "
                                       "Use whatever same as the UCSC genome size file provided; "
                                       "5. Generate md5 sum for each processed file and "
                                       "save to md5_list.txt in the same dir.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_dir",
        type=str,
        required=True,
        help="Path of the ALLC directory."
    )

    parser_req.add_argument(
        "--genome_size_path",
        type=str,
        required=True,
        help="Path of UCSC genome size file, used for check chromosome names. "
             "You should use the same file as the genome you used for mapping "
             "and never change chromosome name by no means."
    )

    parser_opt.add_argument(
        "--compress_level",
        type=int,
        required=False,
        default=6,
        help="Level of ALLC compression level"
    )

    parser_opt.add_argument(
        "--idx",
        type=bool,
        required=False,
        default=True,
        help="Whether add the .idx file for back compatibility."
    )

    parser_opt.add_argument(
        "--remove_additional_chrom",
        type=bool,
        required=False,
        default=False,
        help="Whether remove chroms in ALLC that are not exist in genome size file, "
             "if False and unknown chrom find, will raise error."
    )

    parser_opt.add_argument(
        "--process",
        type=int,
        required=False,
        default=10,
        help="Number of processes to use in parallel."
    )


def mapping_summary_register_subparser(subparser):
    parser = subparser.add_parser('mapping-summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Summary mapping output. Just a convenient function after mapping.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Output directory after mapping."
    )


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'register_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(cemba_data.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True,
                      quiet=False)

    # execute command
    args_vars = vars(args)
    cur_command = args_vars.pop('command')
    # Do real import here:
    if cur_command == 'qsub':
        from .local.qsub import qsub as func
    elif cur_command == 'mapping':
        from .mapping.pipeline import pipeline as func
    elif cur_command == 'default-mapping-config':
        from .mapping.pipeline import print_default_configuration as func
    elif cur_command == 'mapping-qsub':
        from .local.mc.prepare_allc import batch_pipeline as func
    elif cur_command == 'generate-dataset':
        from .local.mc.prepare_dataset import generate_dataset as func
    elif cur_command == 'map-to-region':
        from .tools.allc import map_to_region as func
    elif cur_command == 'assemble-dataset':
        from .local.mc.prepare_dataset import assemble_dataset as func
    elif cur_command == 'allc-to-bigwig':
        from .tools.allc import allc_to_bigwig as func
    elif cur_command == 'merge-allc':
        from .tools.allc_utilities import merge_allc_files as func
    elif cur_command == 'allc-profile':
        from .tools.allc import get_allc_profile as func
    elif cur_command == 'simulate-long-reads-coverage':
        from .tools.simulation import simulate_long_reads_coverage as func
    elif cur_command == 'simulate-allc':
        from .tools.simulation import simulate_allc as func
    elif cur_command == 'allc-extract':
        from .tools.allc import extract_context_allc as func
    elif cur_command == 'allc-standardize':
        from .tools.allc import batch_standardize_allc as func
    elif cur_command == 'mapping-summary':
        from .mapping.pipeline import summary_pipeline_stat as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return


if __name__ == '__main__':
    main()
