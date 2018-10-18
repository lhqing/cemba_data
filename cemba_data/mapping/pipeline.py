import pandas as pd
import glob
import pathlib
import configparser
import os
import logging
import argparse
import json
from datetime import datetime
from .fastq import demultiplex, fastq_qc
from .bismark import bismark
from .allc import call_methylated_sites
from .bam import bam_qc


def _cur_time_string(format='[%Y-%m-%d]%H:%M:%S'):
    return datetime.now().strftime(format)


def _get_configuration(config_path=None):
    ref_path_config = configparser.ConfigParser()
    if config_path is None:
        print('config path not provided, use default config')
        ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    else:
        ref_path_config.read(config_path)
    return ref_path_config


def get_fastq_dataframe(file_path, name_pattern, out_dir=None):
    """
    Generate fastq_dataframe for pipeline input.
    :param file_path: Accept 1. path pattern, 2. path list, 3. path of one file contain all the paths.
    :param name_pattern: Name pattern of file paths, format example: "part1_part2_part3_..._partn",
    should at least contain uid, lane and read_type in the parts. file name are matched to the name patten,
    both should have same length after split by "_".
    :parm out_dir: dir path for fastq dataframe
    :return: fastq_dataframe, required cols: uid, lane, read-type, fastq-path,
    and other optional cols contain metadata.
    """
    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(pathlib.Path(p).absolute()) for p in glob.glob(file_path)]
    elif isinstance(file_path, list):
        pass
    else:
        with open(file_path) as f:
            file_path = [l.strip() for l in f]

    name_pattern = name_pattern.lower().split('_')
    if 'uid' not in name_pattern:
        raise ValueError(f'Name pattern has no "uid"')
    if 'lane' not in name_pattern:
        raise ValueError('Name pattern has no "lane"')
    if 'read-type' not in name_pattern:
        raise ValueError('Name pattern has no "read-type"')

    valid_list = []
    fastq_data = []
    for fastq in file_path:
        fastq_name_list = fastq.split('/')[-1].split('.')[0].split('_')
        if len(fastq_name_list) != len(name_pattern):
            print('Broken name:', fastq)
            continue
        valid_list.append(fastq)
        name_series = pd.Series({p: v for p, v in zip(name_pattern, fastq_name_list)})
        fastq_data.append(name_series)
    fastq_df = pd.DataFrame(fastq_data)
    if '*' in fastq_df.columns:
        del fastq_df['*']
    fastq_df['fastq-path'] = valid_list

    if out_dir is not None:
        fastq_df.to_csv(pathlib.Path(out_dir) / 'fastq_dataframe.tsv.gz',
                        index=None, sep='\t', compression='gzip')

    return fastq_df


def _validate_fastq_dataframe(fastq_dataframe):
    if isinstance(fastq_dataframe, str):
        fastq_dataframe = pd.read_table(fastq_dataframe, index_col=None)
    for required in ['uid', 'lane', 'read-type', 'fastq-path']:
        if required not in fastq_dataframe.columns:
            raise ValueError(required, 'not in fastq dataframe columns')
    return fastq_dataframe


def pipeline(fastq_dataframe, out_dir, config_path=None):
    # get config
    config = _get_configuration(config_path)

    # validate fastq dataframe
    fastq_dataframe = _validate_fastq_dataframe(fastq_dataframe)

    # setup out_dir
    out_dir = pathlib.Path(out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    stat_dir = out_dir / 'stats'
    stat_dir.mkdir()

    # setup logger
    logger = logging.getLogger()
    log_path = stat_dir / "mapping.log"
    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fmt = "%(message)s"
    date_fmt = "[%Y-%m-%d]%H:%M:%S"
    formatter = logging.Formatter(fmt, date_fmt)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # fastq demultiplex
    logger.info('Demultiplex fastq file.')
    demultiplex_df = demultiplex(fastq_dataframe, out_dir, config)
    demultiplex_df.to_csv(stat_dir/'demultiplex_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    logger.info('Trim fastq file and merge lanes.')
    # fastq qc
    fastq_final_df = fastq_qc(demultiplex_df, out_dir, config)
    if fastq_final_df.shape[0] == 0:
        print('no sample remained after fastq qc step')
        return
    else:
        fastq_final_df.to_csv(stat_dir/'fastq_trim_result.tsv.gz',
                              sep='\t', compression='gzip', index=None)

    logger.info('Use bismark and bowtie2 to do mapping.')
    # bismark
    bismark_df = bismark(fastq_final_df, out_dir, config)
    if bismark_df.shape[0] == 0:
        print('no sample remained after bismark step')
        return
    else:
        bismark_df.to_csv(stat_dir / 'bismark_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    logger.info('Deduplicate and filter bam files.')
    # bam
    bam_df = bam_qc(bismark_df, out_dir, config)
    bam_df.to_csv(stat_dir / 'bam_process_result.tsv.gz',
                  sep='\t', compression='gzip', index=None)

    logger.info('Calculate mC sites.')
    # allc
    allc_df = call_methylated_sites(bam_df, out_dir, config)
    allc_df.to_csv(stat_dir / 'allc_total_result.tsv.gz',
                   sep='\t', compression='gzip', index=None)
    return 0


def batch_pipeline(out_dir,
                   fastq_path=None,
                   name_pattern=None,
                   fastq_dataframe_path=None,
                   config_path=None):
    out_dir = pathlib.Path(out_dir).absolute()
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    if config_path is None:
        config_path = os.path.dirname(__file__) + '/mapping_config.ini'
    with open(out_dir / 'mapping_config.ini', 'w') as wf, open(config_path, 'r') as f:
        for line in f:
            wf.write(line)
    config_path = out_dir / 'mapping_config.ini'

    if fastq_path is None and fastq_dataframe_path is None:
        raise ValueError('Both fastq_path and fastq_dataframe_path are None.')
    elif fastq_dataframe_path is not None:
        fastq_dataframe = _validate_fastq_dataframe(fastq_dataframe_path)
    elif fastq_path is not None and name_pattern is None:
        raise ValueError('fastq_path is provided, but name_pattern is None.')
    else:
        fastq_dataframe = get_fastq_dataframe(fastq_path, name_pattern, out_dir=out_dir)

    cmd_list = []
    for uid, sub_df in fastq_dataframe.groupby('uid'):
        uid_out_dir = out_dir / uid
        uid_out_dir.mkdir()
        uid_fastq_dataframe_path = uid_out_dir/'fastq_list.tsv.gz'
        sub_df.to_csv(uid_fastq_dataframe_path,
                      sep='\t', compression='gzip', index=None)
        uid_cmd = f'yap mapping --fastq_dataframe {uid_fastq_dataframe_path} ' \
                  f'--out_dir {uid_out_dir} --config_path {config_path}'
        command_dict = {
            'command': uid_cmd,
            '-pe smp': 20,  # cpu for each command
            # TODO set cup number more cleaver
        }
        cmd_list.append(command_dict)

    with open(out_dir / 'command.json', 'w') as f:
        json.dump(cmd_list, f)


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
        help="Pipeline configuration (.ini) file path"
    )
    return

