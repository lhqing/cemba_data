import logging
import pathlib
import subprocess

import pandas as pd

from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def bismark_bam_qc(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    tmp_dir = output_dir / 'tmp'
    tmp_dir.mkdir(exist_ok=True)

    if isinstance(config, str):
        config = get_configuration(config)
    bismark_records = pd.read_csv(output_dir / 'bismark_mapping.records.csv',
                                  index_col=['uid', 'index_name', 'read_type'],
                                  squeeze=True)
    mapq_threshold = config['bamFilter']['mapq_threshold']

    # process bam
    records = []
    command_list = []
    for (uid, index_name, read_type), bismark_bam_path in bismark_records.iteritems():
        # file path
        filter_bam = bismark_bam_path[:-3] + 'filter.bam'
        sort_bam = bismark_bam_path[:-3] + 'sort.bam'
        final_bam = bismark_bam_path[:-3] + 'final.bam'
        dedup_matrix = bismark_bam_path[:-3] + 'dedup.matrix.txt'
        # command
        filter_cmd = f'samtools view -b -h -q {mapq_threshold} -o {filter_bam} {bismark_bam_path}'
        sort_cmd = f'samtools sort -o {sort_bam} --threads 2 {filter_bam}'
        dedup_cmd = f'picard -Xms4g -Xmx4g MarkDuplicates ' \
                    f'I={sort_bam} O={final_bam} M={dedup_matrix} REMOVE_DUPLICATES=true TMP_DIR={tmp_dir}'
        cleaning_cmd = f'rm -f {bismark_bam_path} {sort_bam} {filter_bam}'
        command = ' && '.join([filter_cmd, sort_cmd, dedup_cmd, cleaning_cmd])
        records.append([uid, index_name, read_type, final_bam])
        command_list.append(command)

    with open(output_dir / 'bismark_bam_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'bam_path'])
    record_df.to_csv(output_dir / 'bismark_bam_qc.records.csv', index=None)
    return record_df, command_list


def summarize_bismark_bam_qc(output_dir):
    bam_dir = pathlib.Path(output_dir)
    output_path = bam_dir / 'bismark_bam_qc.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    bismark_stat_list = list(bam_dir.glob('*.dedup.matrix.txt'))
    for path in bismark_stat_list:
        try:
            report_series = pd.read_csv(path, sep='\t', comment='#').T[0]
        except pd.errors.EmptyDataError:
            # means the bam file is empty
            subprocess.run(['rm', '-f', path])
            continue

        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        report_series['uid'] = uid
        report_series['index_name'] = index_name
        report_series['read_type'] = read_type
        records.append(report_series)
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path, index=None)
    return str(output_path)
