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

log = logging.getLogger()

DESCRIPTION = """
yap (yet another pipeline) is a toolkit for single cell methylation sequencing analysis
Author: Hanqing Liu, hanliu@salk.edu

This toolkit contain functions for 3-stage analysis:
    Stage 1 Preprocessing: Mapping FASTQ, generate single cell ALLC
    Stage 2 Cell Level Analysis: Prepare MCDS dataset for cell based analysis
    Stage 3 Cluster Level Analysis: merge ALLC and ALLC related functions
    Other functions: qsub submitter for SGE QSUB; simulation functions

STAGE 1
    - Mapping pipeline:
    mapping - Actual mapping function
    default-mapping-config - Print the default mapping config
    mapping-qsub - Qsub wrapper for mapping that run with "yap qsub" 
    mapping-summary - Summary mapping output directory after mapping
    
    - Stand alone STAGE 1 functions:
    bam-to-allc - Take 1 after QC and position sorted BAM file, generate and index 1 ALLC file.

STAGE 2
    map-to-region - Map ALLC file into a region BED file
    assemble-dataset - Assemble all cell-region BED file into MCDS
    generate-dataset - Wrapper for map-to-region and assemble-dataset that run with "yap qsub"

STAGE 3
    merge-allc - Merge single cell ALLC files into cluster ALLC.
    allc-profile - Generate summary statistics for a ALLC file.
    allc-to-bigwig - ALLC to BIGWIG.
    allc-extract - Extract ALLC file information. 
    cluster-merge - Wrapper for above stage 3 functions that run with "yap qsub"
    
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

If this looks good, send coffee to...
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
        help="One or space-separated paths of the command file. Accept 2 different command file format: "
             "1. each line in the file is a full command, will be submitted as a single job. "
             "Qsub parameters can only be specified by qsub_global_parms in this way."
             "2. JSON format, a list of dict, where each dict is a command (required key \"command\") "
             "and other optional key specify qsub parameters per job."
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

    parser_opt.add_argument(
        "--wait_until",
        type=str,
        nargs='+',
        required=False,
        help="If provided with a space separate job id(s), "
             "this job will wait until those job finish first."
    )

    parser_opt.add_argument(
        "--force_redo",
        type=bool,
        required=False,
        default=False,
        help="By default, the finished job (which has an accompanying .json file with the command.sh file) "
             "will not be submitted again, but if you want to force resubmit all jobs, set this to true. "
    )

    parser_opt.add_argument(
        "--qsub_global_parms",
        type=str,
        required=False,
        default='',
        help="Other global qsub parameters you want to pass to each job's qsub script. "
             "This will cover command.json if set repeatedly."
             "These additional parameters should form in one ';' separated string like this: "
             "'-q=queue_name;-cwd;-wd=/path/to/working/dir;-pe smp=20;-l h_vmem=5G'"
    )
    return


def pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--fastq_dataframe",
        type=str,
        required=True,
        help="Path of fastq dataframe, can be generate with yap fastq_dataframe"
    )

    parser_req.add_argument(
        "--output_dir",
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

    parser_opt.add_argument(
        "--demultiplex_only",
        dest='demultiplex_only',
        action='store_true',
        help="(Not for mapping) Only demultiplex 8-cell FASTQ into single cell FASTQ by AD index and then quit"
    )
    parser.set_defaults(demultiplex_only=False)

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


def print_plate_info_register_subparser(subparser):
    parser = subparser.add_parser('default-plate-info',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default plate info template.")
    return


def make_sample_sheet_register_subparser(subparser):
    parser = subparser.add_parser('make-sample-sheet',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default plate info template.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--plate_info_paths",
        type=str,
        required=True,
        nargs='+',
        help="Space separated paths of plate infos, at least one file should be provided. "
             "If multiple files provided, will check barcode compatibility."
    )

    parser_req.add_argument(
        "--output_prefix",
        type=str,
        required=True,
        help="Output prefix, will generate 2 sample sheets, 1 for miseq, 1 for novaseq"
    )

    parser_opt.add_argument(
        "--header_path",
        type=str,
        help="Path to the sample sheet header that contains sequencer info. Will use default if not provided."
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
        "--output_dir",
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
        "--output_dir",
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

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="output_dir path that contains all region count bed files, same as map-to-region step"
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
        help="Path of the input ALLC file."
    )

    parser_req.add_argument(
        "--out_prefix",
        type=str,
        required=True,
        help="Path prefix of the output ALLC file."
    )

    parser_req.add_argument(
        "--mc_contexts",
        type=str,
        required=True,
        nargs='+',
        help="Space separated mC contexts to extract."
    )

    parser_opt.add_argument(
        "--strandness",
        type=str,
        required=False,
        default='both',
        help="What to do with strand information, possible values are: "
             "1. both: save +/- strand together in one file without any modification; "
             "2. split: save +/- strand into two separate files, with suffix contain Watson (+) and Crick (-); "
             "3. merge: This will only merge the count on adjacent CpG in +/- strands, only work for CpG like context. "
             "For non-CG context, its the same as both."
    )

    parser_opt.add_argument(
        "--output_format",
        type=str,
        required=False,
        default='allc',
        help="Output format of extracted information, possible values are: "
             "1. allc: keep the allc format; "
             "2. bed5: 5-column bed format, chrom, pos, pos, mc, cov; "
             "3. bg-cov: bedgraph format, chrom, pos, pos, cov; "
             "4. bg-rate: bedgraph format, chrom, pos, pos, mc/cov."
    )

    parser_opt.add_argument(
        "--region",
        type=str,
        required=False,
        default=None,
        help="Only extract records from certain genome region(s) via tabix, "
             "multiple region can be provided in tabix form."
    )


def standardize_allc_register_subparser(subparser):
    parser = subparser.add_parser('allc-standardize',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Standardize one ALLC file by following these rules: "
                                       "1. No header in the ALLC file; "
                                       "2. Chromosome names in ALLC must be same as those "
                                       "in the chrom_size_path file, including 'chr'; "
                                       "3. Output file will be bgzipped with .tbi index; "
                                       "4. Remove additional chromosome (remove_additional_chrom=True) or "
                                       "raise KeyError if unknown chromosome found (default)")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="Path of 1 ALLC file"
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path of UCSC chrom size file"
    )

    parser_opt.add_argument(
        "--compress_level",
        type=int,
        required=False,
        default=5,
        help="Level of bgzip compression"
    )

    parser_opt.add_argument(
        "--idx",
        type=bool,
        required=False,
        default=False,
        help="Whether add the .idx file for back compatibility."
    )

    parser_opt.add_argument(
        "--remove_additional_chrom",
        type=bool,
        required=False,
        default=False,
        help="Whether to remove rows with unknown chromosome instead of raising KeyError"
    )


def mapping_summary_register_subparser(subparser):
    parser = subparser.add_parser('mapping-summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Aggregate mapping summaries.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path of mapping output dir."
    )


def cluster_merge_register_subparser(subparser):
    parser = subparser.add_parser('cluster-merge',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Summary mapping output. Just a convenient function after mapping.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--cluster_table_path",
        type=str,
        required=True,
        help="Path to the cluster table. The first column must be cell id. The other columns "
             "are different levels of cluster assignments. From left to right, sub to major."
    )

    parser_req.add_argument(
        "--allc_path_file",
        type=str,
        required=True,
        help="Path to the ALLC path table. The first column must be cell id. The second column "
             "is the ALLC file path for each cell."
    )

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory, must not exist."
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path to UCSC chrom size file, used for guide merge and bigwig generation."
    )

    parser_opt.add_argument(
        "--bin_length",
        type=int,
        required=False,
        default=1000000,
        help="Length of the chrom bin size when do parallel merging. The larger the more MEM usage. "
             "The default value is usually fine."
    )

    parser_opt.add_argument(
        "--bigwig_contexts",
        type=str,
        required=False,
        nargs='+',
        default=['CGN', 'CHN'],
        help="mC contexts used for generate bigwig."
    )

    parser_opt.add_argument(
        "--extract_contexts",
        type=str,
        required=False,
        nargs='+',
        default=['CGN'],
        help="mC contexts used for extract sub ALLC."
    )

    parser_opt.add_argument(
        "--merge_strand",
        type=bool,
        required=False,
        default=True,
        help="Whether to merge adjacent CG site methylation info."
    )

    parser_opt.add_argument(
        "--min_group",
        type=int,
        required=False,
        default=10,
        help="Minimum number of cells for a cell group to be considered."
    )

    parser_opt.add_argument(
        "--merge_allc_cpu",
        type=int,
        required=False,
        default=20,
        help="CPU used in each merge allc job."
    )


def allc_to_region_count_register_subparser(subparser):
    parser = subparser.add_parser('allc-to-region-count',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Calculate mC and cov at regional level. Region accepted in 2 forms: "
                                       "1. BED file, provided by region_bed_paths, "
                                       "containing arbitrary regions and use bedtools map to calculate; "
                                       "2. Fix-size non-overlap genome bins, provided by bin_sizes, "
                                       "this is much faster to calculate than 1. "
                                       "The output is in 6-column bed-like format: "
                                       "chrom start end region_uid mc cov")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="Path to the ALLC file"
    )

    parser_req.add_argument(
        "--out_prefix",
        type=str,
        required=True,
        help="Path to output prefix"
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path to UCSC chrom size file"
    )

    parser_req.add_argument(
        "--mc_contexts",
        type=str,
        nargs='+',
        required=True,
        help="Space separated mC context list to calculate"
    )

    parser_opt.add_argument(
        "--region_bed_paths",
        type=str,
        nargs='+',
        required=False,
        help="Space separated path list to BED files."
    )

    parser_opt.add_argument(
        "--region_bed_names",
        type=str,
        nargs='+',
        required=False,
        help="Matched name for each BED file provided in region_bed_paths."
    )

    parser_opt.add_argument(
        "--bin_sizes",
        type=int,
        nargs='+',
        required=False,
        help="Space separated genome size bins to calculate."
    )

    parser_opt.add_argument(
        "--max_cov_cutoff",
        type=int,
        required=False,
        default=9999,
        help="Max cov filter for a single site in ALLC"
    )

    parser_opt.add_argument(
        "--save_zero_cov",
        type=bool,
        required=False,
        default=True,
        help="Whether to save the regions that have 0 cov, only apply to region count but not the chromosome count"
    )

    parser_opt.add_argument(
        "--remove_tmp",
        type=bool,
        required=False,
        default=True,
        help="Whether to remove the temporary file"
    )


def bam_to_allc_register_subparser(subparser):
    parser = subparser.add_parser('bam-to-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Take 1 after QC and position sorted BAM file, generate and index 1 ALLC file.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--bam_path",
        type=str,
        required=True,
        help="Path to input position sorted BAM file"
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Path to output ALLC file"
    )

    parser_req.add_argument(
        "--reference_fasta",
        type=str,
        required=True,
        help="Path to 1 genome reference FASTA file (the one used for mapping), "
             "use samtools fadix to build .fai index first. Do not compress that file."
    )

    parser_opt.add_argument(
        "--num_upstr_bases",
        type=int,
        required=False,
        default=0,
        help="Number of upstream base(s) of the C base to include in ALLC context column, "
             "usually use 0 for normal BS-seq, 1 if doing NOMe-seq."
    )

    parser_opt.add_argument(
        "--num_downstr_bases",
        type=int,
        required=False,
        default=2,
        help="Number of downstream base(s) of the C base to include in ALLC context column, "
             "usually use 2 for both BS-seq and NOMe-seq."
    )

    parser_opt.add_argument(
        "--min_mapq",
        type=int,
        required=False,
        default=0,
        help="Minimum MAPQ for a read being considered, samtools mpileup parameter, see samtools documentation. "
             "Redundant if you did that in BAM QC."
    )

    parser_opt.add_argument(
        "--min_base_quality",
        type=int,
        required=False,
        default=1,
        help="Minimum base quality for a base being considered, samtools mpileup parameter, "
             "see samtools documentation. Redundant if you did that in BAM QC."
    )

    parser_opt.add_argument(
        "--compress_level",
        type=int,
        required=False,
        default=5,
        help="Compress level for the output ALLC file, pass to gzip or bgzip"
    )

    parser_opt.add_argument(
        "--idx",
        type=bool,
        required=False,
        default=False,
        help="Whether to generate old methylpy chromosome index, for back compatibility."
    )

    parser_opt.add_argument(
        "--tabix",
        type=bool,
        required=False,
        default=True,
        help="Whether to generate tabix index, bgzip and tabix must be installed for this."
    )

    parser_opt.add_argument(
        "--save_count_df",
        type=bool,
        required=False,
        default=True,
        help="Whether to generate  ALLC context count table."
    )

    parser_opt.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=1,
        help="Number of cores to parallel, NEVER use this if you "
             "run bam-to-allc for things like single cell (generate thousands of ALLCs together)."
    )


def prepare_fastq_dataframe_register_subparser(subparser):
    parser = subparser.add_parser('make-fastq-dataframe',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare fastq_dataframe after prepare_sample_sheet. "
                                       "Only take the fastq files prepared by sample sheet "
                                       "from prepare_sample_sheet function. "
                                       "Otherwise the function may fail or generate wrong file.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--file_path",
        type=str,
        required=True,
        help="Path to input fastq files, either a path pattern contain wildcard "
             "or a path to a file where each row is a fastq path."
    )

    parser_req.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="output path of the fastq dataframe"
    )

    parser_opt.add_argument(
        "--skip_broken_name",
        dest='skip_broken_name',
        action='store_true',
        help='If present, ignore any unrecognized file names in file_path.'
    )
    parser.set_defaults(skip_broken_name=False)


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
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
        from .tools.allc.utilities import map_to_region as func
    elif cur_command == 'assemble-dataset':
        from .local.mc.prepare_dataset import assemble_dataset as func
    elif cur_command == 'allc-to-bigwig':
        from .tools.allc.allc_to_bigwig import allc_to_bigwig as func
    elif cur_command == 'merge-allc':
        from .tools.allc.merge_allc import merge_allc_files as func
    elif cur_command == 'allc-profile':
        from .tools.allc.utilities import profile_allc as func
    elif cur_command == 'simulate-long-reads-coverage':
        from .tools.simulation import simulate_long_reads_coverage as func
    elif cur_command == 'simulate-allc':
        from .tools.simulation import simulate_allc as func
    elif cur_command == 'allc-extract':
        from .tools.allc.extract_allc import extract_allc as func
    elif cur_command == 'allc-standardize':
        from .tools.allc.utilities import standardize_allc as func
    elif cur_command == 'mapping-summary':
        from .mapping.pipeline import summary_pipeline_stat as func
    elif cur_command == 'cluster-merge':
        from .local.mc.prepare_cluster_profile import cluster_merge_pipeline as func
    elif cur_command == 'bam-to-allc':
        from .tools.allc.bam_to_allc import call_methylated_sites as func
    elif cur_command == 'allc-to-region-count':
        from .tools.allc.allc_to_region_count import allc_to_region_count as func
    elif cur_command == 'default-plate-info':
        from .mapping.prepare_sample_sheet import print_plate_info as func
    elif cur_command == 'make-sample-sheet':
        from .mapping.prepare_sample_sheet import make_sample_sheet as func
    elif cur_command == 'make-fastq-dataframe':
        from .local.mc.prepare_allc import get_fastq_dataframe as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"# Executing {cur_command}...")
    func(**args_vars)
    log.info(f"# {cur_command} finished.")
    return


if __name__ == '__main__':
    main()
