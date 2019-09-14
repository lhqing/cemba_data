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
    fastq_dataframe = pd.read_csv(output_dir / 'fastq_dataframe.csv')
    fastq_dataframe = validate_fastq_dataframe(fastq_dataframe)

    random_index_version = config['multiplexIndex']['barcode_version']
    if random_index_version.upper() == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR / 'mapping/files/random_index_v1.fa')
    elif random_index_version.upper() == 'V2':
        # random_index_fasta_path = str(PACKAGE_DIR / 'mapping/files/random_index_v2.fa')
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
        stat_out = output_dir / (f"{uid}-{lane}" + ".demultiplex.stats.txt")
        cmd = f"cutadapt -e 0.01 --no-indels {adapter_parms} " \
              f"-o {r1_out.absolute()} -p {r2_out.absolute()} {r1_in} {r2_in} > {stat_out.absolute()}"
        records.append([uid, lane, str(r1_out), str(r2_out)])
        command_list.append(cmd)
    with open(output_dir / 'demultiplex.command.txt', 'w') as f:
        f.write('\n'.join(command_list))

    # expend records by random index
    # get index names
    index_names = []
    with open(random_index_fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                index_names.append(line.lstrip('>').rstrip())
    expend_records = []
    for record in records:
        for index_name in index_names:
            uid, lane, r1_out, r2_out = record
            expend_records.append([uid, lane, r1_out.format(name=index_name), r2_out.format(name=index_name)])
    record_df = pd.DataFrame(expend_records, columns=['uid', 'lane', 'r1_path', 'r2_path'])
    record_df.to_csv(output_dir / 'demultiplex.records.csv', index=None)
    return record_df, command_list


def _read_cutadapt_result(stat_path):
    """
    Ugly parser of cutadapt output
    TODO: make this nicer, add example output
    """
    with open(stat_path) as f:
        p = re.compile(r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times.")
        series = []
        total_pairs = -1
        for line in f:
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


def summarize_demultiplex(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)
    output_path = output_dir / 'demultiplex.stats.csv'
    if output_path.exists():
        return pd.read_csv(output_path)

    # get index info
    random_index_version = config['multiplexIndex']['barcode_version']
    if random_index_version.upper() == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR / 'mapping/files/random_index_v1.fa')
    elif random_index_version.upper() == 'V2':
        # random_index_fasta_path = str(PACKAGE_DIR / 'mapping/files/random_index_v2.fa')
        # TODO add v2 func and make sure it works
        raise NotImplementedError
    else:
        raise ValueError(f'Unknown version name {random_index_version} in multiplexIndex section of the config file.')
    index_seq_dict = parse_index_fasta(random_index_fasta_path)
    index_name_dict = {v: k for k, v in index_seq_dict.items()}

    # read the demultiplex stats, its per lane, so need to sum up lane together of each uid and index name
    # but R1 R2 is demultiplexed together, so this table don't separate R1 R2
    stat_list = []
    stat_path_list = list(output_dir.glob('*demultiplex.stats.txt'))
    for path in stat_path_list:
        single_df = _read_cutadapt_result(path)
        *uid, suffix = path.name.split('-')
        lane = suffix.split('.')[0]
        uid = '-'.join(uid)
        single_df['uid'] = uid
        single_df['lane'] = lane
        single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
        assert single_df['index_name'].isna().sum() == 0
        stat_list.append(single_df)
    total_demultiplex_stats = pd.concat(stat_list)
    total_demultiplex_stats.to_csv(output_path, index=None)
    subprocess.run(['rm', '-f'] + list(map(str, stat_path_list)))
    return total_demultiplex_stats
