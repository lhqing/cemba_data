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
        default=2500,
        help="Minimum gap distance to be considered as contact"
    )

    parser.add_argument(
        "--keep_split_table",
        dest='keep_split_table',
        action='store_true',
        help='Keep the reads split table before contacts table?'
    )
    parser.set_defaults(keep_split_table=False)


def dss_two_internal_subparser(subparser):
    parser = subparser.add_parser('dss-two',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DSS two-group DMR")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--allc_table_path",
        type=str,
        required=True,
        help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id"
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory"
    )

    parser_req.add_argument(
        "--study_name",
        type=str,
        required=True,
        help="Name of the study"
    )

    parser_req.add_argument(
        "--chrom_sizes_path",
        type=str,
        required=True,
        help="Path to genome chrom size file"
    )

    parser.add_argument(
        "--chroms",
        type=str,
        nargs='+',
        default=None,
        help="Chromosomes to consider."
    )

    parser.add_argument("--p_threshold", type=float, default=0.001,
                        help="FDR threshold to select sig DML")
    parser.add_argument("--min_cg", type=int, required=True, default=1,
                        help="Minimum CpGs for a DMR")
    parser.add_argument("--min_len", type=int, required=True, default=1,
                        help="Minimum length for a DMR")
    parser.add_argument("--sig_ratio", type=float, required=True, default=0.5,
                        help="Minimum ratio of CpGs that are significant in a DMR.")
    parser.add_argument("--delta", type=float, required=True, default=0.1,
                        help="Methylation delta that considered to be informative.")
    parser.add_argument("--cpu", type=int, required=True, default=10,
                        help="Number of CPUs to use")
    parser.add_argument("--chunk_size", type=int, required=False, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument(
        "--not_smoothing",
        dest='smoothing',
        action='store_false',
        help='Do not preform smoothing (will perform smoothing by default).'
    )
    parser.set_defaults(smoothing=True)
    parser.add_argument(
        "--save_dml",
        dest='save_dml',
        action='store_true',
        help='Save all the DML files (will delete DML table by default).'
    )
    parser.set_defaults(save_dml=False)
    return

def dss_multi_internal_subparser(subparser):
    parser = subparser.add_parser('dss-multi',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DSS two-group DMR")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--allc_table_path",
        type=str,
        required=True,
        help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id"
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory"
    )

    parser_req.add_argument(
        "--study_name",
        type=str,
        required=True,
        help="Name of the study"
    )

    parser_req.add_argument(
        "--chrom_sizes_path",
        type=str,
        required=True,
        help="Path to genome chrom size file"
    )

    parser.add_argument(
        "--chroms",
        type=str,
        nargs='+',
        default=None,
        help="Chromosomes to consider."
    )

    parser.add_argument("--p_threshold", type=float, default=0.001,
                        help="FDR threshold to select sig DML")
    parser.add_argument("--min_cg", type=int, required=True, default=1,
                        help="Minimum CpGs for a DMR")
    parser.add_argument("--min_len", type=int, required=True, default=1,
                        help="Minimum length for a DMR")
    parser.add_argument("--sig_ratio", type=float, required=True, default=0.5,
                        help="Minimum ratio of CpGs that are significant in a DMR.")
    parser.add_argument("--cpu", type=int, required=True, default=10,
                        help="Number of CPUs to use")
    parser.add_argument("--chunk_size", type=int, required=False, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument(
        "--no_smooth",
        dest='smoothing',
        action='store_false',
        help='Do not preform smoothing (will perform smoothing by default).'
    )
    parser.set_defaults(smoothing=True)
    parser.add_argument(
        "--save_dml",
        dest='save_dml',
        action='store_true',
        help='Save all the DML files (will delete DML table by default).'
    )
    parser.set_defaults(save_dml=False)
    return


def dmrseq_internal_subparser(subparser):
    parser = subparser.add_parser('dmrseq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Run DMRseq")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--allc_table_path",
        type=str,
        required=True,
        help="Three columns separated by tab: 1) allc path; 2) sample id 3) group id"
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory"
    )

    parser_req.add_argument(
        "--study_name",
        type=str,
        required=True,
        help="Name of the study"
    )

    parser_req.add_argument(
        "--chrom_sizes_path",
        type=str,
        required=True,
        help="Path to genome chrom size file"
    )

    parser.add_argument(
        "--chroms",
        type=str,
        nargs='+',
        default=None,
        help="Chromosomes to consider."
    )

    parser.add_argument("--test_covariate", type=str, default='group',
                        help="Test covariate name in the allc path table")
    parser.add_argument("--match_covariate", type=str, required=False, default=None, nargs='+',
                        help="Matched covariate name in the allc path table (will be adjusted during permutation.")
    parser.add_argument("--adjust_covariate", type=str, required=False, default=None, nargs='+',
                        help="adjust covariate name in the allc path table (will be adjusted in the model.")
    parser.add_argument("--cutoff", type=float, required=False, default=0.1,
                        help="Methylation delta that considered to be informative.")
    parser.add_argument("--min_num_region", type=int, required=False, default=3,
                        help="Minimum number of CpGs that are significant in a DMR, min is 3.")
    parser.add_argument("--bp_span", type=int, required=False, default=1000,
                        help="")
    parser.add_argument("--min_in_span", type=int, required=False, default=30,
                        help="")
    parser.add_argument("--max_gap_smooth", type=int, required=False, default=2500,
                        help="")
    parser.add_argument("--max_gap", type=int, required=False, default=1000,
                        help="")
    parser.add_argument("--max_perms", type=int, required=False, default=10,
                        help="")
    parser.add_argument("--stat", type=str, required=False, default='stat',
                        help="")
    parser.add_argument("--chrs_per_chunk", type=int, required=False, default=1,
                        help="")
    parser.add_argument("--block_size", type=int, required=False, default=5000,
                        help="")
    parser.add_argument("--cpu", type=int, required=False, default=4,
                        help="")
    parser.add_argument("--template_path", type=str, required=False, default='default',
                        help="")
    parser.add_argument("--chunk_size", type=int, required=False, default=50000000,
                        help="chunk size to parallel jobs")
    parser.add_argument(
        "--no_smooth",
        dest='smooth',
        action='store_false',
        help='Do not preform smoothing (will perform smoothing by default).'
    )
    parser.set_defaults(smooth=True)
    parser.add_argument(
        "--block",
        dest='block',
        action='store_true',
        help='Save all the DML files (will delete DML table by default).'
    )
    parser.set_defaults(block=False)
    return


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
    elif cur_command == 'dss-two':
        from .dmr.dss import run_dss_two_group as func
    elif cur_command == 'dss-multi':
        from .dmr.dss import run_dss_multi_group as func
    elif cur_command == 'dmrseq':
        from .dmr.dmrseq import run_dmrseq as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"# Executing {cur_command}...")
    func(**args_vars)
    log.info(f"# {cur_command} finished.")
    return
