import pathlib
import subprocess

import pandas as pd

from .utilities import get_configuration


def generate_allc(input_dir, output_dir, config):
    input_dir = pathlib.Path(input_dir)
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)

    bam_records = pd.read_csv(input_dir / 'final_bam.records.csv',
                              index_col=['uid', 'index_name'],
                              squeeze=True)
    reference_fasta = config['callMethylation']['reference_fasta']
    num_upstr_bases = config['callMethylation']['num_upstr_bases']
    num_downstr_bases = config['callMethylation']['num_downstr_bases']
    compress_level = config['callMethylation']['compress_level']

    records = []
    command_list = []
    for (uid, index_name), bismark_bam_path in bam_records.iteritems():
        # file path
        output_allc_path = output_dir / f'{uid}-{index_name}.allc.tsv.gz'
        # command
        command = f'allcools bam-to-allc ' \
                  f'--bam_path {bismark_bam_path} ' \
                  f'--reference_fasta {reference_fasta} ' \
                  f'--output_path {output_allc_path} ' \
                  f'--cpu 1 ' \
                  f'--num_upstr_bases {num_upstr_bases} ' \
                  f'--num_downstr_bases {num_downstr_bases} ' \
                  f'--compress_level {compress_level} ' \
                  f'--save_count_df'
        records.append([uid, index_name, output_allc_path])
        command_list.append(command)

    with open(output_dir / 'generate_allc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'allc_path'])
    record_df.to_csv(output_dir / 'generate_allc.records.csv', index=None)
    return record_df, command_list


def summarize_generate_allc(output_dir):
    bam_dir = pathlib.Path(output_dir)
    output_path = bam_dir / 'generate_allc.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    allc_stat_list = list(bam_dir.glob('*.count.csv'))
    for path in allc_stat_list:
        try:
            report_df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            # means the bam file is empty
            subprocess.run(['rm', '-f', path])
            continue
        if report_df.shape[0] == 0:
            continue

        *uid, suffix = path.name.split('-')
        uid = '-'.join(uid)
        index_name = suffix.split('.')[0]
        report_df['uid'] = uid
        report_df['index_name'] = index_name
        records.append(report_df)
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.concat(records)
    total_stats_df.to_csv(output_path, index=None)
    return str(output_path)
