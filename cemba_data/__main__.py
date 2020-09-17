"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import logging
import sys

import cemba_data

log = logging.getLogger()

DESCRIPTION = """
YAP (Yet Another Pipeline) is an in-house mapping pipeline for 
snmC-seq, NOMe-seq, snmCT-seq mapping and preprocessing.

See documentation at https://cemba-data.readthedocs.io
"""

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

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the qsub project"
    )

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the qsub project"
    )

    parser.add_argument(
        "--wait_until",
        type=str,
        nargs='+',
        required=False,
        help="If provided with a space separate job id(s), "
             "this job will wait until those job finish first."
    )

    parser.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        default=30,
        help="Total concurrent CPU in qsub list, together with total_mem, "
             "determines how many jobs running in parallel."
    )

    parser.add_argument(
        "--total_mem",
        type=int,
        required=False,
        default=9999,
        help="Total concurrent MEM (GBs) in qsub list, together with total_cpu, "
             "determines how many jobs running in parallel. "
             "Default is 9999, which means only use total_cpu to determine."
    )

    parser.add_argument(
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


def print_plate_info_register_subparser(subparser):
    parser = subparser.add_parser('default-plate-info',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default plate info template.")

    parser.add_argument(
        "--barcode_version", '-v', '-V',
        type=str,
        default='V2',
        choices=['V1', 'V2'],
        help="Use V1 template for 8-random-index library version, "
             "Use V2 template for 384-random-index library version."
    )
    return


def make_sample_sheet_register_subparser(subparser):
    parser = subparser.add_parser('make-sample-sheet',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default plate info template.")

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--plate_info_path",
        "-p",
        type=str,
        required=True,
        help="Path of the plate information file."
    )

    parser_req.add_argument(
        "--output_prefix",
        "-o",
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


def demultiplex_register_subparser(subparser):
    parser = subparser.add_parser('demultiplex',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Demultiplex bcl2fastq results.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--fastq_pattern",
        "-fq",
        type=str,
        required=True,
        help="FASTQ files with wildcard to match all bcl2fastq results, pattern with wildcard must be quoted."
    )

    parser_req.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=True,
        help="Pipeline output directory, will be created recursively."
    )

    parser_req.add_argument(
        "--config_path",
        "-config",
        type=str,
        required=True,
        help="Path to the mapping config, see 'yap default-mapping-config' about how to generate this file."
    )

    parser_req.add_argument(
        "--cpu",
        '-j',
        type=int,
        required=True,
        help="Number of cores to use. Note that the demultiplex step will only use at most 16 cores."
    )
    return


def print_default_config_register_subparser(subparser):
    parser = subparser.add_parser('default-mapping-config',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Print out default config of mapping pipeline")
    from cemba_data.utilities import MAPPING_MODE_CHOICES
    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        choices=MAPPING_MODE_CHOICES,
        help="Library mode"
    )

    parser.add_argument(
        "--barcode_version",
        '-v',
        '-V',
        type=str,
        required=True,
        choices=['V1', 'V2'],
        help="Barcode version, V1 for 8 random index, V2 for 384 random index"
    )

    parser.add_argument(
        "--bismark_ref",
        type=str,
        required=True,
        help="Path to the bismark reference"
    )

    parser.add_argument(
        "--genome_fasta",
        type=str,
        required=True,
        help="Path to the genome fasta file"
    )

    parser.add_argument(
        "--star_ref",
        type=str,
        required=False,
        help="[mct only] Path to the STAR reference, required if mode is mct"
    )

    parser.add_argument(
        "--gtf",
        type=str,
        required=False,
        help="[mct only] Path to the GTF annotation file, required if mode is mct"
    )

    parser.add_argument(
        "--nome",
        dest='nome',
        action='store_true',
        help='Does this library have NOMe treatment?'
    )
    parser.set_defaults(nome=False)
    return


def prepare_register_subparser(subparser):
    parser = subparser.add_parser('prepare',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare batch jobs.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        '-o',
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser.add_argument(
        "--name",
        type=str,
        default=None,
        help="Name of the project, if None, will use the output_dir name."
    )

    parser.add_argument(
        "--total_jobs",
        "-j",
        type=int,
        default=20,
        help="Total number of jobs run in parallel."
    )

    parser.add_argument(
        "--cores_per_job",
        "-c",
        type=int,
        default=5,
        help="Number of cores used by each job, total number of cores is (cores_per_job * total_jobs)."
    )

    parser.add_argument(
        "--memory_per_core",
        "-m",
        type=str,
        default='5G',
        help="Memory assigned to each core, "
             "the total memory of each job is (cores_per_job * memory_per_core); "
             "the total memory of all jobs is (cores_per_job * memory_per_core * total_jobs)"
    )
    return


def mapping_summary_register_subparser(subparser):
    parser = subparser.add_parser('summary',
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
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    qsub_register_subparser(subparsers)
    print_default_config_register_subparser(subparsers)
    print_plate_info_register_subparser(subparsers)
    make_sample_sheet_register_subparser(subparsers)
    mapping_summary_register_subparser(subparsers)
    demultiplex_register_subparser(subparsers)
    prepare_register_subparser(subparsers)


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
    elif cur_command == 'default-plate-info':
        from .demultiplex import print_plate_info as func
    elif cur_command == 'make-sample-sheet':
        from .demultiplex import make_sample_sheet as func
    elif cur_command == 'demultiplex':
        from .demultiplex import demultiplex_pipeline as func
    elif cur_command == 'default-mapping-config':
        from .mapping import print_default_mapping_config as func
    elif cur_command == 'prepare':
        from .mapping import prepare_run as func
    elif cur_command == 'summary':
        from .mapping import final_summary as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    func(**args_vars)
    return


if __name__ == '__main__':
    main()
