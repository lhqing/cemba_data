"""
Input: fasta_dataframe, i7 i5 demultiplexed fastq files by bcl2fastq
Process: remove snmC-seq2 random index by cutadapt to get single cell raw fastq files
Output: Lane merged single cell fastq files

Output file name pattern

r1_out = pathlib.Path(output_dir) / (f"{uid}-{lane}" + "-{name}-R1.fq.gz")
r2_out = pathlib.Path(output_dir) / (f"{uid}-{lane}" + "-{name}-R2.fq.gz")
"""

import locale
import logging
import pathlib
import re
import shlex
import subprocess

import pandas as pd

import cemba_data
from .fastq_dataframe import validate_fastq_dataframe
from .utilities import get_configuration, parse_index_fasta

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def demultiplex(output_dir: str, config: str):
    """
    make a command dataframe that contains all cutadapt demultiplex command
    """
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)

    # read and validate fastq_dataframe
    fastq_dataframe = validate_fastq_dataframe(output_dir / 'fastq_dataframe.csv')

    random_index_version = config['multiplexIndex']['version']
    if random_index_version.upper() == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v1.fa')
    elif random_index_version.upper() == 'V2':
        # random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v2.fa')
        # TODO add v2 func and make sure it works
        raise NotImplementedError
    else:
        raise ValueError(f'Unknown version name {random_index_version} in multiplexIndex section of the config file.')

    # make adapter parms
    adapter_type = '-g'  # The random index is always in R1 5 prime
    adapter_parms = f'{adapter_type} file:{random_index_fasta_path}'

    # standardize read_type
    fastq_dataframe['read_type'] = fastq_dataframe['read_type'].apply(lambda i: 'R2' if '2' in str(i) else 'R1')

    records = []
    command_list = []
    for (uid, lane), sub_df in fastq_dataframe.groupby(['uid', 'lane']):
        tmp_sub_df = sub_df.set_index('read_type')
        r1_in = tmp_sub_df.loc['R1', 'fastq_path']
        r2_in = tmp_sub_df.loc['R2', 'fastq_path']

        # {name} is not format string, its for cutadapt and cutadapt will replace {name} to random index name.
        r1_out = output_dir / (f"{uid}-{lane}" + "-{name}-R1.fq.gz")
        r2_out = output_dir / (f"{uid}-{lane}" + "-{name}-R2.fq.gz")
        cmd = f"cutadapt -e 0.01 --no-indels {adapter_parms} " \
              f"-o {r1_out.absolute()} -p {r2_out.absolute()} {r1_in} {r2_in}"
        records.append([uid, lane, str(r1_out).replace('{name}', '*'), str(r2_out).replace('{name}', '*')])
        command_list.append(cmd)
    with open(output_dir / 'demultiplex.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records, columns=['uid', 'lane', 'r1_path_pattern', 'r2_path_pattern'])
    record_df.to_csv(output_dir / 'demultiplex.records.csv', index=None)
    return record_df, command_list


def _read_cutadapt_result(result):
    """
    Ugly parser of cutadapt output
    TODO: make this nicer, add example output
    """
    p = re.compile(r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times.")
    series = []
    total_pairs = -1
    for line in result.split('\n'):
        if line.startswith('Total read pairs processed'):
            # some weird transform of cutadapt outputs...
            locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
            total_pairs = locale.atoi(line.split(' ')[-1])

        m = p.search(line)
        if m is not None:
            result_dict = {}
            for i in m.group().split('; '):
                k, v = i.split(': ')
                result_dict[k] = v
            result_series = pd.Series(result_dict)
            series.append(result_series)

    total_df = pd.DataFrame(series)
    total_df['Trimmed'] = total_df['Trimmed'].apply(lambda c: c.split(' ')[0]).astype(int)
    total_df['TotalPair'] = total_pairs
    total_df['Ratio'] = total_df['Trimmed'] / total_pairs
    return total_df


def demultiplex_runner(command):
    """
    cutadapt command wrapper, run command and parse results
    """
    try:
        cutadapt_run = subprocess.run(shlex.split(command),
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE,
                                      encoding='utf8',
                                      check=True)
        result_df = _read_cutadapt_result(cutadapt_run.stdout)
    except subprocess.CalledProcessError as e:
        log.error("Pipeline break, FASTQ demultiplex ERROR! command was:")
        log.error(command)
        log.error(e.stdout)
        log.error(e.stderr)
        raise e

    p = re.compile(r'-o (?P<r1_output_path>.+-R1.fq.gz)')
    r1_output_path = pathlib.Path(p.search(command)['r1_output_path'])
    stat_out = r1_output_path.parent / ('-'.join(r1_output_path.name.split('-')[:4]) + '-demultiplex.stat.csv')
    result_df.to_csv(stat_out)
    return


def summarize_demultiplex(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)

    # get index info
    random_index_version = config['multiplexIndex']['version']
    if random_index_version.upper() == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v1.fa')
    elif random_index_version.upper() == 'V2':
        # random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v2.fa')
        # TODO add v2 func and make sure it works
        raise NotImplementedError
    else:
        raise ValueError(f'Unknown version name {random_index_version} in multiplexIndex section of the config file.')
    index_seq_dict = parse_index_fasta(random_index_fasta_path)
    index_name_dict = {v: k for k, v in index_seq_dict.items()}

    # read the demultiplex stats, its per lane, so need to sum up lane together of each uid and index name
    # but R1 R2 is demultiplexed together, so this table don't separate R1 R2
    stat_list = []
    stat_path_list = list(output_dir.glob('*demultiplex.stat.csv'))
    for path in stat_path_list:
        single_df = pd.read_csv(path, index_col=0)
        *uid, lane, _ = path.name.split('-')
        uid = '-'.join(uid)
        single_df['uid'] = uid
        single_df['lane'] = lane
        single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
        assert single_df['index_name'].isna().sum() == 0
        stat_list.append(single_df)
    total_demultiplex_stats = pd.concat(stat_list)
    total_demultiplex_stats.to_csv(output_dir / 'demultiplex.stat.total.csv', index=None)
    subprocess.run(['rm', '-f'] + list(map(str, stat_path_list)))
    return total_demultiplex_stats
