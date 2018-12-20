import argparse
import sys
import inspect
import cemba_data
import logging
from .description import description, epilog

log = logging.getLogger()


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
        help="Path of the command dict file"
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        help="Total CPU in qsub list"
    )

    parser_opt.add_argument(
        "--submission_gap",
        type=int,
        required=False,
        help="Submission Gap in qsub list"
    )

    parser_opt.add_argument(
        "--qstat_gap",
        type=int,
        required=False,
        help="Qstat check gap in qsub list"
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


def batch_map_to_region_register_subparser(subparser):
    parser = subparser.add_parser('map-to-region-qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")
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

    return


def assemble_dataset_register_subparser(subparser):
    parser = subparser.add_parser('assemble-dataset',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")

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
        "--dataset_name",
        type=str,
        required=True,
        help="Name of the dataset"
    )

    parser_opt.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=10,
        help="Number of cores to use."
    )


def merge_allc_register_subparser(subparser):
    parser = subparser.add_parser('merge-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Just a wrapper of methylpy merge_allc_files, "
                                       "without doing methylpy's index. "
                                       "But use bgzip and tabix instead")

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

    parser_opt.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=1,
        help="Number of CPUs for merge ALLC, parallel on chrom bins level."
    )

    parser_opt.add_argument(
        "--index",
        type=bool,
        required=False,
        default=False,
        help="methylpy default index, not doing it by default."
    )

    parser_opt.add_argument(
        "--get_mcg",
        type=bool,
        required=False,
        default=True,
        help="Generate a CG only, strand merged ALLC file for CG DMR calling."
    )

    parser_opt.add_argument(
        "--cg_pattern",
        type=str,
        required=False,
        default='CGN',
        help="mCG context pattern for --get_mcg."
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


def main():
    parser = argparse.ArgumentParser(description=description,
                                     epilog=epilog)
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
    elif cur_command == 'map-to-region-qsub':
        from .local.mc.prepare_dataset import batch_map_to_region as func
    elif cur_command == 'map-to-region':
        from .tools.allc import map_to_region as func
    elif cur_command == 'assemble_dataset':
        from .local.mc.prepare_dataset import assemble_dataset as func
    elif cur_command == 'allc-to-bigwig':
        from .tools.allc import allc_to_bigwig as func
    elif cur_command == 'merge-allc':
        from .tools.allc import merge_allc as func
    else:
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return


if __name__ == '__main__':
    main()
