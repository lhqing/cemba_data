import pandas as pd
import glob
import pathlib
import configparser
import os


def _get_configuration():
    ref_path_config = configparser.ConfigParser()
    ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    return ref_path_config


def _get_fastq_dataframe(file_path, name_pattern):
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


def pipeline(fastq_dataframe, out_dir):
    # validate fastq dataframe
    fastq_dataframe = validate_fastq_dataframe(fastq_dataframe)

    # setup out_dir
    out_dir = pathlib.Path(out_dir).absolute()
    table_dir = out_dir / 'table'
    fastq_dir = out_dir / 'fastq'
    bam_dir = out_dir / 'bam'
    allc_dir = out_dir / 'allc'
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
        table_dir.mkdir()
        fastq_dir.mkdir()
        bam_dir.mkdir()
        allc_dir.mkdir()

    #





