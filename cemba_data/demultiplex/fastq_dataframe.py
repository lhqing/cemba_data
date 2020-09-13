"""
Generate raw FASTQ dataframe based on fixed name pattern
name pattern is based on samplesheet generated in plateinfo_and_samplesheet.py
"""

import glob
import logging
import pathlib

import pandas as pd

# logger
log = logging.getLogger()


def _parse_v1_fastq_path(path):
    """
    UID pattern of V1 {sample_id_prefix}-{plate1}-{plate2}-{plate_pos}
    FASTQ name pattern of V1:
    {sample_id_prefix}-{plate1}-{plate2}-{plate_pos}_{internal_info}_{lane}_{read_type}_{internal_info}.fastq.gz
    """
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
        raise ValueError(f'Found unknown name pattern in path {path}')
    name_dict = dict(plate1=plate1,
                     plate2=plate2,
                     plate_pos=plate_pos,
                     lane=lane,
                     read_type=read_type,
                     fastq_path=path,
                     uid=f'{plate1}-{plate2}-{plate_pos}')
    name_series = pd.Series(name_dict)
    return name_series


def _parse_v2_fastq_path(path):
    """
    UID pattern of V2 {sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}
    FASTQ name pattern of V1:
    {sample_id_prefix}-{plate}-{multiplex_group}-{barcode_name}_{internal_info}_{lane}_{read_type}_{internal_info}.fastq.gz
    """
    path = pathlib.Path(path)
    try:
        *_, plate, multiplex_group, multi_field = path.name.split('-')
        primer_name, _, lane, read_type, _ = multi_field.split('_')
        try:
            assert primer_name[0] in 'ABCDEFGHIJKLMNOP'
            assert int(primer_name[1:]) in list(range(1, 25))
            assert int(multiplex_group) in list(range(1, 7))
            assert lane in {'L001', 'L002', 'L003', 'L004'}
            assert read_type in {'R1', 'R2'}
        except AssertionError:
            raise ValueError
    except ValueError:
        raise ValueError(f'Found unknown name pattern in path {path}')
    name_dict = dict(plate=plate,
                     multiplex_group=multiplex_group,
                     primer_name=primer_name,
                     lane=lane,
                     read_type=read_type,
                     fastq_path=path,
                     uid=f'{plate}-{multiplex_group}-{primer_name}')
    name_series = pd.Series(name_dict)
    return name_series


def make_fastq_dataframe(file_path, barcode_version, output_path=None):
    """
    Generate fastq_dataframe for pipeline input.

    Parameters
    ----------
    file_path
        Accept 1. path pattern contain wildcard, 2. path list, 3. path of one file contain all the paths.
    barcode_version
        Only accept two options: 1) V1 for 8 random index; 2) V2 for 384 random index.
    output_path
        output path of the fastq dataframe
    Returns
    -------
        fastq_dataframe for pipeline input.
    """
    barcode_version = barcode_version.upper()
    if barcode_version == 'V1':
        parser = _parse_v1_fastq_path
    elif barcode_version == 'V2':
        parser = _parse_v2_fastq_path
    else:
        raise ValueError(f'Primer Version can only be V1 or V2, got {barcode_version}.')

    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(pathlib.Path(p).absolute()) for p in glob.glob(file_path)]
    elif isinstance(file_path, list):
        pass
    else:
        with open(file_path) as f:
            file_path = [line.strip() for line in f]
    log.info(f'{len(file_path)} FASTQ file paths in input')

    fastq_data = []
    for path in file_path:
        name_series = parser(path)
        fastq_data.append(name_series)
    fastq_df = pd.DataFrame(fastq_data)
    log.info(f'{fastq_df.shape[0]} valid fastq names.')
    if fastq_df.shape[0] == 0:
        log.info('No fastq name remained, check if the name pattern is correct.')
        return None

    # make sure UID is unique
    for _, df in fastq_df.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError(f'UID column is not unique.')
    if output_path is not None:
        fastq_df.to_csv(output_path, index=False)
    return fastq_df
