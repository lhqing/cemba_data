import glob
import logging
import pathlib

import pandas as pd

# logger
log = logging.getLogger()


def get_fastq_dataframe(file_path, output_path=None, skip_broken_name=False):
    """
    Generate fastq_dataframe for pipeline input.

    Parameters
    ----------
    file_path
        Accept 1. path pattern contain wildcard, 2. path list, 3. path of one file contain all the paths.
    output_path
        output path of the fastq dataframe
    skip_broken_name
        If true, ignore any unrecognized file names in file_path
    Returns
    -------
        fastq_dataframe for pipeline input.
    """
    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(pathlib.Path(p).absolute()) for p in glob.glob(file_path)]
    elif isinstance(file_path, list):
        pass
    else:
        with open(file_path) as f:
            file_path = [l.strip() for l in f]
    log.info(f'{len(file_path)} fastq file paths in input')

    fastq_data = []
    broken_names = []
    for path in file_path:
        path = pathlib.Path(path)
        try:
            *_, plate1, plate2, multi_field = path.name.split('-')
            plate_pos, _, lane, read_type, _ = multi_field.split('_')
            try:
                assert plate_pos[0] in 'ABCDEFGH'
                assert int(plate_pos[1:]) in list(range(1, 13))
                assert lane in {'L001', 'L002', 'L003', 'L004'}
                assert read_type in {'R1', 'R2'}
                assert plate1 != plate2
            except AssertionError:
                raise ValueError
        except ValueError:
            broken_names.append(path)
            if skip_broken_name:
                continue
            raise ValueError(f'Found unknown name pattern in path {path}')
        name_dict = dict(plate1=plate1,
                         plate2=plate2,
                         plate_pos=plate_pos,
                         lane=lane,
                         read_type=read_type,
                         fastq_path=path,
                         uid=f'{plate1}-{plate2}-{plate_pos}')
        name_series = pd.Series(name_dict)
        fastq_data.append(name_series)

    fastq_df = pd.DataFrame(fastq_data)
    log.info(f'{len(broken_names)} broken names.')
    log.info(f'{fastq_df.shape[0]} valid fastq names.')
    if fastq_df.shape[0] == 0:
        log.info('No fastq name remained, check if the name pattern is correct.')
        return None

    # make sure UID is unique
    for _, df in fastq_df.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError(f'UID column is not unique.')
    if output_path is not None:
        if not output_path.endswith('tsv.gz'):
            output_path += 'tsv.gz'
        fastq_df.to_csv(output_path, index=None, sep='\t', compression='gzip')
        return
    else:
        return fastq_df


def validate_fastq_dataframe(fastq_dataframe):
    """
    Check if fastq_dataframe is
    1. have required columns
    2. uid is unique
    """
    if isinstance(fastq_dataframe, str):
        fastq_dataframe = pd.read_csv(fastq_dataframe, index_col=None, sep='\t')

    for required in ['uid', 'lane', 'read_type', 'fastq_path']:
        if required not in fastq_dataframe.columns:
            raise ValueError(f'column {required} not in fastq dataframe columns, '
                             f'remember that the 4 required columns of fastq dataframe are: '
                             f'uid, lane, read_type, fastq_path. '
                             f'The orders do not matter, but the names need to be exact.')

    for _, df in fastq_dataframe.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError('uid column are not unique for each lane and read-type combination.')

    # modify fastq dataframe column names
    fastq_dataframe.columns = [column.replace('-', '_') for column in fastq_dataframe.columns]

    # modify fastq columns, because '_' is used in file name and we split by '_'
    # I know this is stupid...
    fastq_dataframe['uid'] = fastq_dataframe['uid'].apply(lambda i: i.replace('_', '-'))
    fastq_dataframe['lane'] = fastq_dataframe['lane'].apply(lambda i: i.replace('_', '-'))
    fastq_dataframe['read_type'] = fastq_dataframe['read_type'].apply(lambda i: i.replace('_', '-'))
    return fastq_dataframe
