import logging
import pathlib
import subprocess

import pandas as pd

import cemba_data
from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def fastq_qc(output_dir, config):
    # TODO add pair end fastq_qc for STAR pair end mapping
    """
    reads level QC and trimming. R1 R2 separately.
    """
    output_dir = pathlib.Path(output_dir)
    merge_records_df = pd.read_csv(output_dir / 'merge_lane.records.csv').set_index(['uid', 'index_name'])
    if isinstance(config, str):
        config = get_configuration(config)

    # parameters
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

    # determine whether proceed based on number of trimmed reads
    use_pairs = set()
    try:
        total_demultiplex_stats = pd.read_csv(output_dir / 'demultiplex.stats.csv')
        for (uid, index_name), sub_df in total_demultiplex_stats.groupby(['uid', 'index_name']):
            sample_demultiplex_total = sub_df['Trimmed'].sum()
            if sample_demultiplex_total < 1:
                print(f'In  uid {uid}: index {index_name} is empty.')
                continue
            if sample_demultiplex_total < total_reads_threshold:
                print(f'In  uid {uid}: index {index_name} skipped '
                      f'due to too less reads: {sample_demultiplex_total}')
                continue
            use_pairs.add((uid, index_name))
    except FileNotFoundError:
        # stats file not found, could due to dry run mode.
        # In this case, the command will be generate for every pair of cell
        for (uid, index_name), sub_df in merge_records_df.groupby(['uid', 'index_name']):
            use_pairs.add((uid, index_name))

    records = []
    command_list = []
    for (uid, index_name, read_type), sub_df in merge_records_df.groupby(['uid', 'index_name', 'read_type']):
        if (uid, index_name) not in use_pairs:
            continue
        fastq_in_path = sub_df['fastq_path'][0]

        if read_type.upper() == 'R1':
            _left_cut = r1_left_cut
            _right_cut = r1_right_cut
            _adapter = r1_adapter
        elif read_type.upper() == 'R2':
            _left_cut = r2_left_cut
            _right_cut = r2_right_cut
            _adapter = r2_adapter
        else:
            raise ValueError(f'Unknown read type {read_type}.')

        output_path = f'{output_dir}/{uid}-{index_name}-{read_type}.trimmed.fq.gz'
        cmd = f'cutadapt -j {cutadapt_cores} --report=minimal -O {overlap} ' \
              f'-q {quality_threshold} -u {_left_cut} ' \
              f'-u -{_right_cut} -m {length_threshold} ' \
              f'-a {_adapter} -o {output_path} {fastq_in_path} > {output_path}.stats.txt'
        records.append([uid, index_name, read_type, output_path])
        command_list.append(cmd)

    with open(output_dir / 'fastq_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'fastq_path'])
    record_df.to_csv(output_dir / 'fastq_qc.records.csv', index=None)
    return record_df, command_list


def summarize_fastq_qc(output_dir):
    fastq_dir = pathlib.Path(output_dir)
    output_path = fastq_dir / 'fastq_qc.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    fastq_stat_list = list(fastq_dir.glob('*.fastq_qc.stats.csv'))
    for path in fastq_stat_list:
        fastq_qc_stat = pd.read_csv(path, index_col=0, squeeze=True)
        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        fastq_qc_stat['uid'] = uid
        fastq_qc_stat['index_name'] = index_name
        fastq_qc_stat['read_type'] = read_type
        records.append(fastq_qc_stat)
        subprocess.run(['rm', '-f', path], check=True)
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path)
    return str(output_path)
