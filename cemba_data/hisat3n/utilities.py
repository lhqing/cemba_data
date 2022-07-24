import pathlib
import re
import yaml
import pandas as pd

from ..utilities import get_configuration


def _read_yaml_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def _read_ini_config(config_path):
    return get_configuration(config_path)


def read_mapping_config(cwd: str = '.'):
    tried = []
    yaml_path = None
    for name in ['config', 'mapping_config']:
        for config_dir in [cwd, f'{cwd}/..']:
            for suffix in ['yaml', 'yml']:
                path = f'{config_dir}/{name}.{suffix}'
                tried.append(path)
                if pathlib.Path(path).exists():
                    yaml_path = path
    default_path = f'~/mapping_config.yaml'
    if pathlib.Path(default_path).exists():
        yaml_path = default_path

    ini_path = None
    for name in ['config', 'mapping_config']:
        for config_dir in [cwd, f'{cwd}/..']:
            path = f'{config_dir}/{name}.ini'
            tried.append(path)
            if pathlib.Path(path).exists():
                ini_path = path

    if yaml_path is not None:
        config = _read_yaml_config(yaml_path)
    elif ini_path is not None:
        config = _read_ini_config(ini_path)
    else:
        config = {}
    return config


def validate_cwd_fastq_paths(cwd: str = '.'):
    """
    Validate fastq paths in the fastq subdirectory of cwd.
    Parameters
    ----------
    cwd :
        Path of the current working directory.

    Returns
    -------
    fastq_table : pandas.DataFrame
    """
    # get all fastq file paths
    fastq_paths = [p
                   for p in pathlib.Path(f'{cwd}/fastq/').glob('*.[fq.gz][fastq.gz]')
                   if 'trim' not in p.name]

    # parse cell id and match fastq pairs
    fastq_pattern = re.compile(r'(?P<cell_id>.+)(-|_)(?P<read_type>(R1|R2|r1|r2)).(fastq|fq)(.gz)*')
    fastq_records = {}
    for p in fastq_paths:
        match = fastq_pattern.match(p.name)
        if match is None:
            # print(f'WARNING: {p} has FASTQ file path suffix, but do not match '
            #       f'expected file name pattern {fastq_pattern}')
            pass
        else:
            cell_id = match.group('cell_id')
            read_type = match.group('read_type')
            fastq_records[cell_id, read_type.upper()] = str(p)

    if len(fastq_records) == 0:
        raise ValueError('No fastq files found in fastq folder, '
                         'or no fastq files match expected file name pattern')

    fastq_table = pd.Series(fastq_records).unstack()
    if 'R1' not in fastq_table.columns or 'R2' not in fastq_table.columns:
        raise ValueError('No R1 or R2 fastq files found')
    fastq_table = fastq_table[['R1', 'R2']].copy()

    # raise error if fastq file not paired
    missing_file = fastq_table.isna().sum(axis=1) > 0
    if missing_file.sum() > 0:
        for cell in missing_file[missing_file].index:
            print(f'{cell} missing R1 or R2 FASTQ file.')
        raise FileNotFoundError(f'FASTQ files in {pathlib.Path(f"{cwd}/fastq/").absolute()} is not all paired.')
    return fastq_table
