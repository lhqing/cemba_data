"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import inspect
import logging
import sys

import cemba_data
from .__main__ import setup_logging

log = logging.getLogger()

DESCRIPTION = """
yap-internal is used for automation, not intend to be used by end user. 
Use yap instead. 
"""

EPILOG = ''


def select_dna_reads_internal_subparser(subparser):
    parser = subparser.add_parser('select-dna-reads',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Select DNA reads, see cemba_data.mapping.mct_bismark_bam_filter")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--input_bam",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_bam",
        type=str,
        required=True
    )

    parser.add_argument(
        "--mc_rate_max_threshold",
        type=float,
        default=0.5,
    )

    parser.add_argument(
        "--cov_min_threshold",
        type=int,
        default=5
    )
    return


def select_rna_reads_internal_subparser(subparser):
    parser = subparser.add_parser('select-rna-reads',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Select RNA reads, see cemba_data.mapping.mct_star_bam_filter")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--input_bam",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--output_bam",
        type=str,
        required=True
    )

    parser.add_argument(
        "--mc_rate_min_threshold",
        type=float,
        default=0.9,
    )

    parser.add_argument(
        "--cov_min_threshold",
        type=int,
        default=5
    )
    return


def featurecount_internal_subparser(subparser):
    parser = subparser.add_parser('featurecount',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Batch FeatureCount wrapper")
    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--bam_table",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--out_prefix",
        type=str,
        required=True
    )

    parser_req.add_argument(
        "--gtf_path",
        type=str,
        required=True
    )

    parser.add_argument(
        "--count_type",
        type=str,
        default='gene',
    )

    parser.add_argument(
        "--id_type",
        type=str,
        default='gene_id',
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=2
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=50
    )
    return


def atac_bulk_pipeline_internal_subparser(subparser):
    parser = subparser.add_parser('atac-bulk-pipeline',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Simple Wrapper. "
                                       "Make cluster level merged fragment bed, bigwig "
                                       "and MACS2 peaks from ATAC clustering results.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--cell_group_path",
        type=str,
        required=True,
        help="first row is header, names of cluster columns will be used as output names; "
             "Column 0 is sample name; "
             "Column 1 is cell barcode; "
             "Column 2, ..., n are level of cell group/clusters;"
    )

    parser_req.add_argument(
        "--output_dir_path",
        type=str,
        required=True,
        help="Output directory, each cluster col will be a sub-dir"
    )

    parser_req.add_argument(
        "--sample_snap_path",
        type=str,
        required=True,
        help="no header; "
             "Column 0 is sample name; "
             "Column 1 is sample SNAP file path;"
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="chromosome size file path"
    )

    parser_req.add_argument(
        "--species",
        type=str,
        required=True,
        help="hs or mm"
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=1,
        help="Number of cpu to parallel"
    )


def mapping_summary_internal_subparser(subparser):
    parser = subparser.add_parser('summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="mapping summary CLI for internal ues")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="output_dir of each UID"
    )


def split_read_internal_subparser(subparser):
    parser = subparser.add_parser('m3c-split-reads',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Split unmapped reads for remap in snm3C-seq data")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--fastq_path",
        type=str,
        required=True,
        help="Input fastq path"
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Output fastq path"
    )

    parser_req.add_argument(
        "--trim_b",
        type=int,
        default=0,
        help="Whether trim on both ends before splitting reads. Default is not trim."
    )

    parser_req.add_argument(
        "--size_l",
        type=int,
        default=40,
        help="Size to slice on the left part"
    )

    parser_req.add_argument(
        "--size_r",
        type=int,
        default=40,
        help="Size to slice on the right part"
    )

    parser_req.add_argument(
        "--size_m",
        type=int,
        default=30,
        help="Minimum size of the middle part"
    )


def generate_contacts_internal_subparser(subparser):
    parser = subparser.add_parser('generate-contacts',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="parse ")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--bam_path",
        type=str,
        required=True,
        help="Input bam path"
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Output contact file path"
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Chrom size path"
    )

    parser_req.add_argument(
        "--min_gap",
        type=int,
        default=1000,
        help="Minimum gap distance to be considered as contact"
    )


def internal_main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'internal_subparser' in name:
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
    if cur_command == 'select-dna-reads':
        from .mapping.mct import select_dna_reads as func
    elif cur_command == 'select-rna-reads':
        from .mapping.mct import select_rna_reads as func
    elif cur_command == 'summary':
        from .mapping.stats import mapping_stats as func
    elif cur_command == 'm3c-split-reads':
        from .mapping.m3c import split_fastq_reads as func
    elif cur_command == 'generate-contacts':
        from .mapping.m3c import generate_contacts as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"# Executing {cur_command}...")
    func(**args_vars)
    log.info(f"# {cur_command} finished.")
    return
