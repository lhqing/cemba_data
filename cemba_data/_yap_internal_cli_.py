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

    parser.add_argument(
        "--remove_input",
        dest='remove_input',
        action='store_true')
    parser.set_defaults(remove_input=False)
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

    parser.add_argument(
        "--remove_input",
        dest='remove_input',
        action='store_true')
    parser.set_defaults(remove_input=False)
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
        from .mapping import select_dna_reads as func
    elif cur_command == 'select-rna-reads':
        from .mapping import select_rna_reads as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"# Executing {cur_command}...")
    func(**args_vars)
    log.info(f"# {cur_command} finished.")
    return
