import pathlib
import subprocess
from collections import defaultdict

import pandas as pd
import pysam
from ALLCools._open import open_bam

from .utilities import get_bam_header_str, get_configuration

METHYLATED_CHAR = 'H'
UNMETHYLATED_CHAR = 'h'


def read_mc_level(bismark_tag):
    m_c = sum([bismark_tag.count(c) for c in METHYLATED_CHAR])
    normal_c = sum([bismark_tag.count(c) for c in UNMETHYLATED_CHAR])
    total_c = m_c + normal_c
    if total_c == 0:
        return 0, 0
    else:
        read_mc_rate = m_c / total_c
        return read_mc_rate, total_c


def select_dna_reads(input_bam,
                     output_bam,
                     mc_rate_max_threshold=0.5,
                     cov_min_threshold=5,
                     remove_input=True):
    bam_header = get_bam_header_str(input_bam)
    read_profile_dict = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as f, open_bam(output_bam, 'w') as out_f:
        out_f.write(bam_header)
        for read in f:
            bismark_tag = read.get_tag('XM')
            mc_rate, cov = read_mc_level(bismark_tag)
            read_profile_dict[(int(100 * mc_rate), cov)] += 1

            # split reads
            if (mc_rate > mc_rate_max_threshold) or (cov < cov_min_threshold):
                continue
            out_f.write(read.tostring() + '\n')
    read_profile = pd.Series(read_profile_dict)
    read_profile.index.name = ['mc_rate', 'cov']
    read_profile.to_csv(str(output_bam) + '.reads_profile.csv', header=True)
    if remove_input:
        subprocess.run(['rm', '-f', input_bam])
    return


def prepare_select_dna_reads(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    if isinstance(config, str):
        config = get_configuration(config)

    bismark_records = pd.read_csv(output_dir / 'bismark_bam_qc.records.csv',
                                  index_col=['uid', 'index_name', 'read_type'],
                                  squeeze=True)
    mc_rate_max_threshold = config['DNAReadsFilter']['mc_rate_max_threshold']
    cov_min_threshold = config['DNAReadsFilter']['cov_min_threshold']
    remove_input = config['DNAReadsFilter']['remove_input']

    # process bam
    records = []
    command_list = []
    for (uid, index_name, read_type), bismark_bam_path in bismark_records.iteritems():
        # file path
        output_bam = bismark_bam_path[:-3] + 'dna_reads.bam'
        # command
        keep_input_str = '--remove_input' if remove_input else ''
        command = f'yap-internal select-dna-reads ' \
                  f'--input_bam {bismark_bam_path} ' \
                  f'--output_bam {output_bam} ' \
                  f'--mc_rate_max_threshold {mc_rate_max_threshold} ' \
                  f'--cov_min_threshold {cov_min_threshold} ' \
                  f'{keep_input_str}'
        records.append([uid, index_name, read_type, output_bam])
        command_list.append(command)

    with open(output_dir / 'select_dna_reads.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'bam_path'])
    record_df.to_csv(output_dir / 'select_dna_reads.records.csv', index=None)
    return record_df, command_list


def summarize_select_dna_reads(output_dir):
    bam_dir = pathlib.Path(output_dir)
    output_path = bam_dir / 'select_dna_reads.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    select_dna_reads_stat_list = list(bam_dir.glob('*.reads_profile.csv'))
    for path in select_dna_reads_stat_list:
        try:
            report_df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            # means the bam file is empty
            subprocess.run(['rm', '-f', path])
            continue
        if report_df.shape[0] == 0:
            subprocess.run(['rm', '-f', path])
            continue

        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        report_df['uid'] = uid
        report_df['index_name'] = index_name
        report_df['read_type'] = read_type
        records.append(report_df)
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.concat(records)
    total_stats_df.to_csv(output_path, index=None)
    return str(output_path)
