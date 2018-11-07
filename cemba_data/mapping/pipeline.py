import pandas as pd
import pathlib
import configparser
import os
import argparse
from datetime import datetime
from .fastq import demultiplex, fastq_qc
from .bismark import bismark
from .allc import call_methylated_sites
from .bam import bam_qc


def _cur_time_string(formats='[%Y-%m-%d]%H:%M:%S'):
    return datetime.now().strftime(formats)


def _get_configuration(config_path=None):
    ref_path_config = configparser.ConfigParser()
    if config_path is None:
        print('config path not provided, use default config')
        ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    else:
        ref_path_config.read(config_path)
    return ref_path_config


def print_default_configuration(out_path=None):
    with open(os.path.dirname(__file__) + '/mapping_config.ini') as f:
        configs = f.readlines()
    if out_path is not None:
        with open(out_path, 'w') as f:
            f.writelines(configs)
    else:
        for line in configs:
            print(line, end='')
    return


def validate_fastq_dataframe(fastq_dataframe):
    if isinstance(fastq_dataframe, str):
        fastq_dataframe = pd.read_table(fastq_dataframe, index_col=None)

    for required in ['uid', 'lane', 'read_type', 'fastq_path']:
        if required not in fastq_dataframe.columns:
            raise ValueError(required, 'not in fastq dataframe columns')

    for _, df in fastq_dataframe.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError(f'uid column are not unique for each lane and read-type combination.')

    # modify fastq dataframe column names
    fastq_dataframe.columns = [column.replace('-', '_') for column in fastq_dataframe.columns]

    # modify fastq columns, because '_' is used in file name and we split by '_'
    fastq_dataframe['uid'] = fastq_dataframe['uid'].apply(lambda i: i.replace('_', '-'))
    fastq_dataframe['lane'] = fastq_dataframe['lane'].apply(lambda i: i.replace('_', '-'))
    fastq_dataframe['read_type'] = fastq_dataframe['read_type'].apply(lambda i: i.replace('_', '-'))

    return fastq_dataframe


def pipeline(fastq_dataframe, out_dir, config_path=None):
    # get config
    config = _get_configuration(config_path)

    # get and validate fastq dataframe
    fastq_dataframe = validate_fastq_dataframe(fastq_dataframe)

    # setup out_dir
    out_dir = pathlib.Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    stat_dir = out_dir / 'stats'
    stat_dir.mkdir()

    # setup logger
    # TODO add log module, replace all print

    # fastq demultiplex
    print('Demultiplex fastq file.')
    demultiplex_df = demultiplex(fastq_dataframe, out_dir, config)
    demultiplex_df.to_csv(stat_dir/'demultiplex_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    print('Trim fastq file and merge lanes.')
    # fastq qc
    fastq_final_df = fastq_qc(demultiplex_df, out_dir, config)
    if fastq_final_df.shape[0] == 0:
        print('no sample remained after fastq qc step')
        return
    else:
        fastq_final_df.to_csv(stat_dir/'fastq_trim_result.tsv.gz',
                              sep='\t', compression='gzip', index=None)

    print('Use bismark and bowtie2 to do mapping.')
    # bismark
    bismark_df = bismark(fastq_final_df, out_dir, config)
    if bismark_df.shape[0] == 0:
        print('no sample remained after bismark step')
        return
    else:
        bismark_df.to_csv(stat_dir / 'bismark_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    print('Deduplicate and filter bam files.')
    # bam
    bam_df = bam_qc(bismark_df, out_dir, config)
    bam_df.to_csv(stat_dir / 'bam_process_result.tsv.gz',
                  sep='\t', compression='gzip', index=None)

    print('Calculate mC sites.')
    # allc
    allc_df = call_methylated_sites(bam_df, out_dir, config)
    allc_df.to_csv(stat_dir / 'allc_total_result.tsv.gz',
                   sep='\t', compression='gzip', index=None)

    print('Mapping finished.')
    return 0


def pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")
    parser.set_defaults(func=pipeline)

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
    parser.set_defaults(func=print_default_configuration)

    parser_opt = parser.add_argument_group("Optional inputs")

    parser_opt.add_argument(
        "--out_path",
        type=str,
        required=False,
        help="Path to save the config file"
    )
    return
