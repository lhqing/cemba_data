"""
Input: fastq dataframe, output fastq files from bcl2fastq
Processes:
    - remove adapter
    - demultiplex fastq by AD index
    - trim reads by base quality, uniform cut and filter by length
    - merge lane
Output: fastq_final_result dataframe
"""

import pathlib
import pandas as pd
import subprocess
import re
import locale
import functools
import multiprocessing
import shlex
import logging
from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def _make_command_dataframe(fastq_dataframe, out_dir, config):
    """
    make a command dataframe that contains all cutadapt demultiplex command
    """
    multiplex_index_dict = config['multiplexIndex']
    overlap = int(config['demultiplex']['overlap'])
    anchor = bool(config['demultiplex']['anchor'])
    adapter_pos = int(config['demultiplex']['adapter_pos'])

    # make adapter parms
    adapter_type = '-g' if int(adapter_pos) == 5 else '-a'
    if anchor:
        # X stands for anchored adapter trim
        multiplex_index_dict = {k: 'X' + v for k, v in multiplex_index_dict.items()}
    adapter_parms = ' '.join([f'{adapter_type} {k}={v}' for k, v in multiplex_index_dict.items()])

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
        cmd = f"cutadapt {adapter_parms} -O {overlap} -o {r1_out.absolute()} -p {r2_out.absolute()} {r1_in} {r2_in}"
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
        pipeline universal output_dir
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


def fastq_qc(demultiplex_result, out_dir, config):
    """
    reads level QC and trimming. R1 R2 separately and merge Lane together.

    Parameters
    ----------
    demultiplex_result
        dataframe from demultiplex step
    out_dir
        pipeline universal output_dir
    config
        pipeline universal config
    Returns
    -------
    fastq_final_result
        id columns: uid, index_name, read_type
    """

    if isinstance(config, str):
        config = get_configuration(config)

    cutadapt_cores = int(config['fastqTrim']['cutadapt_cores'])

    r1_adapter = config['fastqTrim']['r1_adapter']
    r2_adapter = config['fastqTrim']['r1_adapter']
    length_threshold = config['fastqTrim']['length_threshold']
    quality_threshold = config['fastqTrim']['quality_threshold']
    r1_left_cut = config['fastqTrim']['r1_left_cut']
    r1_right_cut = config['fastqTrim']['r1_right_cut']
    r2_left_cut = config['fastqTrim']['r2_left_cut']
    r2_right_cut = config['fastqTrim']['r2_right_cut']
    overlap = config['fastqTrim']['overlap']
    total_reads_threshold = int(config['fastqTrim']['total_reads_threshold'])

    results = []
    for (uid, index_name), sub_df in demultiplex_result.groupby(['uid', 'index_name']):
        sample_demultiplex_total = sub_df['Trimmed'].sum()
        if sample_demultiplex_total < 1:
            log.info(f'In  uid {uid}: index {index_name} is empty.')
            continue
        if sample_demultiplex_total < total_reads_threshold:
            log.info(f'In  uid {uid}: index {index_name} skipped '
                     f'due to too less reads: {sample_demultiplex_total}')
            continue
        # merge R1
        r1_path_pattern = f'{out_dir}/{uid}_*_{index_name}_R1.fq.gz'
        r1_merge_cmd = f'pigz -cd -p {cutadapt_cores} {r1_path_pattern} | pigz'
        r1_merge_result = subprocess.run(r1_merge_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         shell=True, check=True)
        r1_raw_out = f'{out_dir}/{uid}_{index_name}_R1.raw.fq.gz'
        with open(r1_raw_out, 'wb') as f:
            f.write(r1_merge_result.stdout)

        # process R1
        r1_out = f'{out_dir}/{uid}_{index_name}_R1.trimed.fq.gz'
        r1_cmd = f'cutadapt -j {cutadapt_cores} --report=minimal -O {overlap} ' \
                 f'-q {quality_threshold} -u {r1_left_cut} ' \
                 f'-u -{r1_right_cut} -m {length_threshold} ' \
                 f'-a {r1_adapter} -o {r1_out} {r1_raw_out}'
        try:
            r1_result = subprocess.run(r1_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       encoding='utf8', shell=True, check=True)

        except subprocess.CalledProcessError as e:
            log.error("Pipeline break, FASTQ R1 trim ERROR!")
            log.error(e.stdout)
            log.error(e.stderr)
            raise e

        # get R1 result stat
        lines = []
        for line in r1_result.stdout.split('\n'):
            ll = line.split('\t')
            if len(ll) > 1:
                lines.append(ll)
        s = pd.Series({name: number for name, number in zip(*lines)})
        s['uid'] = uid
        s['index_name'] = index_name
        s['read_type'] = 'R1'
        results.append(s)

        # merge R2
        r2_path_pattern = f'{out_dir}/{uid}_*_{index_name}_R2.fq.gz'
        r2_merge_cmd = f'pigz -cd -p {cutadapt_cores} {r2_path_pattern} | pigz'
        r2_merge_result = subprocess.run(r2_merge_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                         shell=True, check=True)
        r2_raw_out = f'{out_dir}/{uid}_{index_name}_R2.raw.fq.gz'
        with open(r2_raw_out, 'wb') as f:
            f.write(r2_merge_result.stdout)

        # process R2
        r2_out = f'{out_dir}/{uid}_{index_name}_R2.trimed.fq.gz'
        r2_cmd = f'cutadapt -j {cutadapt_cores} --report=minimal -O {overlap} ' \
                 f'-q {quality_threshold} -u {r2_left_cut} ' \
                 f'-u -{r2_right_cut} -m {length_threshold} ' \
                 f'-a {r2_adapter} -o {r2_out} {r2_raw_out}'
        try:
            r2_result = subprocess.run(r2_cmd, stdout=subprocess.PIPE,
                                       encoding='utf8', shell=True, check=True)
        except subprocess.CalledProcessError as e:
            log.error("Pipeline break, FASTQ R2 trim ERROR!")
            log.error(e.stdout)
            log.error(e.stderr)
            raise e

        # get R2 result stat
        lines = []
        for line in r2_result.stdout.split('\n'):
            ll = line.split('\t')
            if len(ll) > 1:
                lines.append(ll)
        s = pd.Series({name: number for name, number in zip(*lines)})
        s['uid'] = uid
        s['index_name'] = index_name
        s['read_type'] = 'R2'
        results.append(s)

    fastq_final_result = pd.DataFrame(results)
    if len(results) == 0:
        # all sample skipped
        return fastq_final_result
    fastq_final_result['out_reads_rate'] = \
        fastq_final_result['out_reads'].astype(int) / fastq_final_result['in_reads'].astype(int)
    fastq_final_result['out_bp_rate'] = \
        fastq_final_result['out_reads'].astype(int) / fastq_final_result['in_reads'].astype(int)

    # clean up
    for (uid, index_name), sub_df in demultiplex_result.groupby(['uid', 'index_name']):
        r_path_pattern = f'{out_dir}/{uid}_*_{index_name}_R*.fq.gz'
        r_rm_cmd = f'ionice -c 2 -n 0 rm -f {r_path_pattern}'
        subprocess.run(r_rm_cmd, shell=True)
    for uid in demultiplex_result['uid'].unique():
        # remove unknown reads
        r_path_pattern = f'{out_dir}/{uid}_L*_unknown_R*.fq.gz'
        r_rm_cmd = f'ionice -c 2 -n 0 rm -f {r_path_pattern}'
        subprocess.run(r_rm_cmd, shell=True)

    return fastq_final_result
