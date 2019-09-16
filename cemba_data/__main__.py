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


def print_default_config_register_subparser(subparser):
    parser = subparser.add_parser('default-mapping-config',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default config of mapping pipeline")
    from cemba_data.mapping.pipeline import MAPPING_MODE_CHOICES
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=MAPPING_MODE_CHOICES,
        help="Library mode"
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


def pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Batch mapping pipeline from bcl2fastq multiplexed FASTQ file to ALLC file. ")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--input_fastq_pattern",
        type=str,
        required=True,
        help="File path pattern with wildcards of demultiplexed fastq files"
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

    parser.add_argument(
        "--fastq_dataframe_path",
        type=str,
        default=None,
        help="FASTQ dataframe path, this is optional, "
             "if the library is demultiplexed using SampleSheet generated by yap make-sample-sheet, "
             "FASTQ dataframe will be generated automatically."
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=5,
        help="Number of cores to use in different steps, the effect depending on --mode: "
             "if --mode command_only, --cpu has no effect; "
             "if --mode qsub, --cpu is the sum of the slots active qsub job occupy; "
             "if --mode local, --cpu is the number of process to use in current machine"
    )

    parser.add_argument(
        "--mode",
        type=str,
        default='command_only',
        choices=['command_only', 'qsub', 'local'],
        help="Run mode, if command_only, will only generate command files but not execute; "
             "if qsub, will execute with SGE qsub system; "
             "if local, will execute with current system, only use this for debugging."
    )
    return


def mapping_summary_register_subparser(subparser):
    parser = subparser.add_parser('mapping-summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping summary after the pipeline finished.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )


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
        from cemba_data.qsub import qsub as func
    elif cur_command == 'default-mapping-config':
        from .mapping.pipeline import print_default_configuration as func
    elif cur_command == 'default-plate-info':
        from .mapping import print_plate_info as func
    elif cur_command == 'make-sample-sheet':
        from .mapping import make_sample_sheet as func
    elif cur_command == 'mapping':
        from .mapping.pipeline import pipeline as func
    elif cur_command == 'mapping-summary':
        from .mapping.summary import mapping_summary as func
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
