import pandas as pd
import glob
import pathlib
import configparser
import os
from .fastq import _demultiplex, _fastq_qc
from .bismark import _bismark
from .allc import _call_methylated_sites
from .bam import _bam_qc


def _get_configuration(config_path=None):
    ref_path_config = configparser.ConfigParser()
    if config_path is None:
        print('config path not provided, use default config')
        ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    else:
        ref_path_config.read(config_path)
    return ref_path_config


def get_fastq_dataframe(file_path, name_pattern):
    """
    Generate fastq_dataframe for pipeline input.
    :param file_path: Accept 1. path pattern, 2. path list, 3. path of one file contain all the paths.
    :param name_pattern: Name pattern of file paths, format example: "part1_part2_part3_..._partn",
    should at least contain uid, lane and read_type in the parts. file name are matched to the name patten,
    both should have same length after split by "_".
    :return: fastq_dataframe, required cols: uid, lane, read-type, fastq-path,
    and other optional cols contain metadata.
    """
    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(p.absolute()) for p in glob.glob(file_path)]
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
    return fastq_df


def validate_fastq_dataframe(fastq_dataframe):
    if isinstance(fastq_dataframe, str):
        fastq_dataframe = pd.read_table(fastq_dataframe, index_col=None)
    for required in ['uid', 'lane', 'read-type', 'fastq-path']:
        if required not in fastq_dataframe.columns:
            raise ValueError(required, 'not in fastq dataframe columns')

    # TODO add other validations

    return fastq_dataframe


def pipeline(fastq_dataframe, out_dir, config_path=None):
    # get config
    config = _get_configuration(config_path)

    # validate fastq dataframe
    fastq_dataframe = validate_fastq_dataframe(fastq_dataframe)

    # setup out_dir
    out_dir = pathlib.Path(out_dir).absolute()
    stat_dir = out_dir / 'stats'
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
        stat_dir.mkdir()

    # fastq demultiplex
    demultiplex_df = _demultiplex(fastq_dataframe, out_dir, config)
    demultiplex_df.to_csv(stat_dir/'demultiplex_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    # fastq qc
    fastq_final_df = _fastq_qc(demultiplex_df, out_dir, config)
    fastq_final_df.to_csv(stat_dir/'fastq_trim_result.tsv.gz',
                          sep='\t', compression='gzip', index=None)

    # bismark
    bismark_df = _bismark(fastq_final_df, out_dir, config)
    bismark_df.to_csv(stat_dir / 'bismark_result.tsv.gz',
                      sep='\t', compression='gzip', index=None)

    # bam
    bam_df = _bam_qc(bismark_df, out_dir, config)
    bam_df.to_csv(stat_dir / 'bam_process_result.tsv.gz',
                  sep='\t', compression='gzip', index=None)

    # allc
    allc_df = _call_methylated_sites(bam_df, out_dir, config)
    allc_df.to_csv(stat_dir / 'allc_total_result.tsv.gz',
                   sep='\t', compression='gzip', index=None)
