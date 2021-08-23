"""
When adding new function:
1. add a func_register_subparser function to register the subparser
2. add a condition in main func about this new func name, import the real func as func in main
"""

import argparse
import logging
import sys

import cemba_data
from cemba_data import __version__

log = logging.getLogger()

DESCRIPTION = """
YAP (Yet Another Pipeline) is a mapping pipeline for multiple
snmC-seq based single-cell sequencing technologies.

See documentation at https://hq-1.gitbook.io/mc/
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


def sbatch_register_subparser(subparser):
    parser = subparser.add_parser('sbatch',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="General sbatch helper, "
                                       "need to prepare a command dict file where key is a command, "
                                       "value is a list of target files.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the qsub project"
    )

    parser_req.add_argument(
        "--command_file_path",
        type=str,
        required=True,
        help="Each row is a command."
    )

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the qsub project"
    )

    parser_req.add_argument(
        "--time_str",
        type=str,
        required=True,
        help="Time estimate (upper limit) of Sbatch Jobs."
    )

    parser_req.add_argument(
        "--queue",
        type=str,
        default='skx-normal',
        help="Queue partition of stampede2. "
             "See https://portal.tacc.utexas.edu/user-guides/stampede2#running-queues for more information."
    )

    parser_req.add_argument(
        "--email",
        type=str,
        default=None,
        help="If you want to get notice about job status, put your email here."
    )

    parser_req.add_argument(
        "--email_type",
        type=str,
        default='fail',
        help="Type of status you want to get notice, default is fail, means only get notice when job failed."
    )

    parser.add_argument(
        "--max_jobs",
        type=int,
        required=False,
        default=None,
        help="Max number of jobs for the same user (not for the same sbatch command). "
             "If not provided, will determine automatically based on stampede2 limits."
    )

    parser.add_argument(
        "--dry_run",
        required=False,
        action='store_true',
        help="Prepare scripts for sbatch manual submission or debug."
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
        help="Number of cores to use. Note that the demultiplex step will only use at most 32 cores, "
             "the merge lane step will use the number of cores you provided."
    )
    return


def update_snakemake_register_subparser(subparser):
    parser = subparser.add_parser('update-snakemake',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Update Snakemake related files when mapping_config.ini "
                                       "OR file directory changed.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=True,
        help="Pipeline output directory."
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
        "--chrom_size_path",
        type=str,
        required=False,
        help="[m3c only] Path to the chrom size file, "
             "only chromosomes occur in this file will be considered in generating chromatin contacts, "
             "required if mode is mct"
    )

    parser.add_argument(
        "--nome",
        dest='nome',
        action='store_true',
        help='Does this library have NOMe treatment?'
    )
    parser.set_defaults(nome=False)
    return


def start_from_cell_fastq_register_subparser(subparser):
    parser = subparser.add_parser('start-from-cell-fastq',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare batch jobs directly from cell-level FASTQ files, "
                                       "no demultiplex step.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        '-o',
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser_req.add_argument(
        "--config_path",
        "-config",
        type=str,
        required=True,
        help="Path to the mapping config, see 'yap default-mapping-config' about how to generate this file."
    )

    parser_req.add_argument(
        "--fastq_pattern",
        "-fq",
        type=str,
        required=True,
        help="Path pattern with wildcard to match all cell-level FASTQ files, pattern with wildcard must be quoted."
    )
    return


def summary_register_subparser(subparser):
    parser = subparser.add_parser('summary',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping summary after the pipeline finished.")

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--output_dir",
        '-o',
        type=str,
        required=True,
        help="Pipeline output directory."
    )

    parser_req.add_argument(
        "--notebook",
        '-nb',
        type=str,
        required=False,
        default=None,
        help="Notebook template for mapping summary, if not provided, will use yap default template."
    )
    return


def mc_bulk_subparser(subparser):
    parser = subparser.add_parser('mc-bulk',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare the snakefile for merging single-cell ALLC files "
                                       "into pseudo-bulk ALLC files and generate BigWig files.")
    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_table",
        "-p",
        type=str,
        required=True,
        help="Path of the allc table. The allc table is a two column tsv file. "
             "The first columns is the absolute ALLC file paths; "
             "the second column is the group name of each file."
    )

    parser_req.add_argument(
        "--output_dir",
        "-o",
        type=str,
        required=True,
        help="Path of the output directory, will be created if not exist."
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path of the chromosome size file path."
    )

    parser_opt.add_argument(
        "--mch_context",
        type=str,
        default='CHN',
        help="mCH contexts for generating the bigwig tracks."
    )

    parser_opt.add_argument(
        "--mcg_context",
        type=str,
        default='CGN',
        help="mCG contexts for generating the bigwig tracks and merge strand."
    )

    parser_opt.add_argument(
        "--bigwig_mch_bin_size",
        type=int,
        default=50,
        help="Bin size used to generate mCH bigwig."
    )

    parser_opt.add_argument(
        "--bigwig_mcg_bin_size",
        type=int,
        default=1,
        help="Bin size used to generate mCG bigwig."
    )

    parser_opt.add_argument(
        "--cpu_per_job",
        type=int,
        default=12,
        help="Number of CPUs to use in individual merge-allc job."
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        default=12,
        help="Number of CPUs to use in total."
    )


def snm3c_imputation_subparser(subparser):
    parser = subparser.add_parser('m3c-impute',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare the snakefile for snm3C-seq contact matrices imputation.")
    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path of the output directory after mapping is finished."
    )

    parser_req.add_argument(
        "--chrom_size_path",
        type=str,
        required=True,
        help="Path of the chromosome sizes file."
    )

    parser_opt.add_argument(
        "--scheduler",
        type=str,
        default='sbatch',
        choices=['sbatch', 'qsub'],
        help="Which scheduler to use, if using sbatch, remember to provide remote chrom_size_path."
    )

    parser_opt.add_argument(
        "--scool_cpu",
        type=int,
        default=5,
        help="Number of CPU to use in generating scool files."
    )

    parser_opt.add_argument(
        "--cpu_per_job",
        type=int,
        default=10,
        help="Number of CPU for each imputation job."
    )

    parser_opt.add_argument(
        "--batch_size",
        type=int,
        default=100,
        help="Number of cells to process in each imputation job."
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        default=100,
        help="Number of total CPUs to use in qsub."
    )

    parser_opt.add_argument(
        "--min_contacts_per_cell",
        type=int,
        default=50000,
        help="Minimum contacts for a cell to be included in the contacts analysis."
    )

    parser_opt.add_argument(
        "--sbatch_time_str",
        type=str,
        default='2:00:00',
        help="Time estimation for sbatch jobs."
    )

    parser_opt.add_argument(
        "--min_cutoff",
        type=str,
        default='1e-5',
        help="Minimum value cutoff of the imputation matrix. "
             "Values smaller than the cutoff will be set to 0 to reduce file size."
    )

    parser_opt.add_argument(
        "--skip_scool_prep",
        dest='skip_scool_prep',
        action='store_true',
        help='Skip scool preparation step if you already have them prepared in the same location '
             'and just want to use this command to change settings and update snakefiles.'
    )
    parser.set_defaults(skip_scool_prep=False)
    return


def snm3c_dataset_subparser(subparser):
    parser = subparser.add_parser('m3c-dataset',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare the snakefile for snm3C-seq contact matrices imputation.")
    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path of the output directory after mapping is finished."
    )

    parser_req.add_argument(
        "--fasta_path",
        type=str,
        required=True,
        help="Path of the chromosome fasta file."
    )

    parser_opt.add_argument(
        "--cpu",
        type=int,
        default=10,
        help="Number of CPU for each job."
    )


def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version", help="Show version number and exit",
                        version=__version__)
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    qsub_register_subparser(subparsers)
    sbatch_register_subparser(subparsers)
    print_default_config_register_subparser(subparsers)
    print_plate_info_register_subparser(subparsers)
    make_sample_sheet_register_subparser(subparsers)
    demultiplex_register_subparser(subparsers)
    update_snakemake_register_subparser(subparsers)
    start_from_cell_fastq_register_subparser(subparsers)
    summary_register_subparser(subparsers)
    mc_bulk_subparser(subparsers)
    snm3c_imputation_subparser(subparsers)
    snm3c_dataset_subparser(subparsers)

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
    elif cur_command == 'sbatch':
        from cemba_data.sbatch import sbatch_submitter as func
    elif cur_command == 'default-plate-info':
        from .demultiplex import print_plate_info as func
    elif cur_command == 'make-sample-sheet':
        from .demultiplex import make_sample_sheet as func
    elif cur_command == 'update-snakemake':
        from .demultiplex import update_snakemake as func
    elif cur_command == 'demultiplex':
        from .demultiplex import demultiplex_pipeline as func
    elif cur_command == 'default-mapping-config':
        from .mapping import print_default_mapping_config as func
    elif cur_command == 'start-from-cell-fastq':
        from .mapping import start_from_cell_fastq as func
    elif cur_command == 'summary':
        from cemba_data.mapping import final_summary as func
    elif cur_command == 'mc-bulk':
        from cemba_data.bulk import prepare_mc_bulk as func
    elif cur_command == 'm3c-impute':
        from cemba_data.snm3C import prepare_impute_dir as func
    elif cur_command == 'm3c-dataset':
        from cemba_data.snm3C import prepare_dataset_commands as func
    else:
        log.debug(f'{cur_command} not Known, check the main function if else part')
        parser.parse_args(["-h"])
        return

    # run the command
    func(**args_vars)
    return


if __name__ == '__main__':
    main()
