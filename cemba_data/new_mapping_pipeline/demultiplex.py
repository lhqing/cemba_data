"""
Input: i7 i5 demultiplexed fastq files by bcl2fastq
Process: remove snmC-seq2 random index by cutadapt to get single cell raw fastq files
Output: Lane merged single cell fastq files
"""

import functools
import locale
import logging
import multiprocessing
import pathlib
import re
import shlex
import subprocess

import pandas as pd

import cemba_data

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def _make_command_dataframe(fastq_dataframe, out_dir, config):
    """
    make a command dataframe that contains all cutadapt demultiplex command
    """
    random_index_version = config['multiplexIndex']['version']
    if random_index_version.upper() == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v1.fa')
    elif random_index_version.upper() == 'V2':
        random_index_fasta_path = str(PACKAGE_DIR / 'new_mapping_pipeline/files/random_index_v2.fa')
        # TODO add v2 func and make sure it works
        raise NotImplementedError
    else:
        raise ValueError(f'Unknown version name {random_index_version} in multiplexIndex section of the config file.')

    overlap = int(config['demultiplex']['overlap'])
    adapter_pos = int(config['demultiplex']['adapter_pos'])

    # make adapter parms
    adapter_type = '-g' if int(adapter_pos) == 5 else '-a'
    adapter_parms = f'{adapter_type} file:{random_index_fasta_path}'

    # make cmd_list
    required_cols = ('uid', 'lane', 'read_type', 'fastq_path')
    for col in required_cols:
        if col not in fastq_dataframe.columns:
            raise ValueError(col, 'not in fastq dataframe')
    # standardize read_type
    fastq_dataframe['read_type'] = fastq_dataframe['read_type'].apply(lambda i: 'R2' if '2' in str(i) else 'R1')

    records = []
    for (uid, lane), sub_df in fastq_dataframe.groupby(['uid', 'lane']):
        tmp_sub_df = sub_df.set_index('read_type')
        r1_in = tmp_sub_df.loc['R1', 'fastq_path']
        r2_in = tmp_sub_df.loc['R2', 'fastq_path']
        r1_out = pathlib.Path(out_dir) / (f"{uid}_{lane}" + "_{name}_R1.fq.gz")
        r2_out = pathlib.Path(out_dir) / (f"{uid}_{lane}" + "_{name}_R2.fq.gz")
        cmd = f"cutadapt -e 0.01 --no-indels {adapter_parms} " \
              f"-O {overlap} -o {r1_out.absolute()} -p {r2_out.absolute()} {r1_in} {r2_in}"
        records.append([uid, lane, cmd])

    cmd_df = pd.DataFrame(records, columns=['uid', 'lane', 'cmd'])
    return cmd_df


def _read_cutadapt_result(result):
    """
    Ugly parser of cutadapt output
    TODO: make this nicer, add example output
    """
    p = re.compile(r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times.")
    series = []
    total_pairs = -1
    for line in result.stdout.split('\n'):
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


def demultiplex(fastq_dataframe, out_dir, config):
    """
    demultiplex of AD index using cutadapt. R1 R2 together, and each lane separately.

    Parameters
    ----------
    fastq_dataframe
        pipeline input fastq_dataframe
    out_dir
        pipeline universal out_dir
    config
        pipeline universal config

    Returns
    -------
    demultiplex result dataframe, id columns: uid, index_name, lane.

    """

    if isinstance(config, str):
        config = get_configuration(config)

    multiplex_index_dict = config['multiplexIndex']
    multiplex_index_map = {v: k for k, v in multiplex_index_dict.items()}
    cmd_df = _make_command_dataframe(fastq_dataframe, out_dir, config)

    cutadapt_run = functools.partial(subprocess.run,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     encoding='utf8',
                                     check=True)
    total_results = []
    for uid, single_uid_df in cmd_df.groupby('uid'):
        pool = multiprocessing.Pool(cmd_df.shape[0])
        results = []
        for i, row in single_uid_df.iterrows():
            result = pool.apply_async(cutadapt_run, (shlex.split(row['cmd']),))
            results.append(result)
        pool.close()
        pool.join()

        for result, lane in zip(results, single_uid_df['lane']):
            try:
                result_df = _read_cutadapt_result(result.get())
            except subprocess.CalledProcessError as e:
                log.error("Pipeline break, FASTQ demultiplex ERROR!")
                log.error(e.stdout)
                log.error(e.stderr)
                raise e

            result_df['lane'] = lane
            result_df['uid'] = uid
            result_df['index_name'] = result_df['Sequence'].apply(lambda ind: multiplex_index_map[ind])
            total_results.append(result_df)
    total_result_df = pd.concat(total_results, ignore_index=True)
    return total_result_df
