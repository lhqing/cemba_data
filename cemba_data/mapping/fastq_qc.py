import logging
import pathlib
import subprocess

import pandas as pd

import cemba_data
from .utilities import get_configuration, parse_index_fasta

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def fastq_qc(output_dir, config):
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
    total_demultiplex_stats.to_csv(output_dir / 'demultiplex.stat.total.csv')
    subprocess.run(['rm', '-f'] + list(map(str, stat_path_list)))

    # determine whether proceed based on number of trimmed reads
    use_pairs = set()
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
              f'-a {_adapter} -o {output_path} {fastq_in_path}'
        records.append([uid, index_name, read_type, output_path])
        command_list.append(cmd)

    with open(output_dir / 'fastq_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'fastq_path'])
    record_df.to_csv(output_dir / 'fastq_qc.records.csv', index=None)
    return


def fastq_qc_runner(command):
    try:
        p = subprocess.run(command,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           encoding='utf8',
                           shell=True,
                           check=True)
        # get R1 result stat
        output_path = command.split(' ')[-2]
        lines = []
        for line in p.stdout.split('\n'):
            ll = line.split('\t')
            if len(ll) > 1:
                lines.append(ll)
        s = pd.Series({name: number for name, number in zip(*lines)})
        s.to_csv(output_path + '.fastq_qc.stats.csv', header=True)
    except subprocess.CalledProcessError as e:
        log.error("Got error in fastq_qc, command was:")
        log.error(command)
        log.error(e.stdout)
        log.error(e.stderr)
        raise e


def summarize_fastq_qc(fastq_dir):
    fastq_dir = pathlib.Path(fastq_dir)
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
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path)
    return str(output_path)
