import logging
import pathlib

import pandas as pd

from .bismark_mapping import summarize_bismark_mapping
from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def bismark_bam_qc(output_dir, config):
    """
    """
    output_dir = pathlib.Path(output_dir)
    summarize_bismark_mapping(bam_dir=output_dir)
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
        sort_bam = bismark_bam_path[:-3] + 'sort.bam'
        dedup_bam = bismark_bam_path[:-3] + 'dedup.bam'
        dedup_matrix = bismark_bam_path[:-3] + 'dedup.matrix.txt'
        filter_bam = bismark_bam_path[:-3] + 'filter.bam'
        # command
        sort_cmd = f'samtools sort -o {sort_bam} --threads 2 {bismark_bam_path}'
        dedup_cmd = f'picard MarkDuplicates I={sort_bam} O={dedup_bam} M={dedup_matrix} REMOVE_DUPLICATES=true'
        filter_cmd = f'samtools view -b -h -q {mapq_threshold} -o {filter_bam} {dedup_bam}'
        cleaning_cmd = f'rm -f {bismark_bam_path} {sort_bam} {dedup_bam}'
        command = ' && '.join([sort_cmd, dedup_cmd, filter_cmd, cleaning_cmd])
        records.append([uid, index_name, read_type, filter_bam])
        command_list.append(command)

    with open(output_dir / 'bismark_bam_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'bam_path'])
    record_df.to_csv(output_dir / 'bismark_bam_qc.records.csv', index=None)
    return


def star_bam_qc(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    if isinstance(config, str):
        config = get_configuration(config)
    star_records = pd.read_csv(output_dir / 'star_mapping.records.csv',
                               index_col=['uid', 'index_name'],
                               squeeze=True)

    mapq_threshold = config['bamFilter']['mapq_threshold']

    # process bam
    records = []
    command_list = []
    for (uid, index_name), star_bam_path in star_records.iteritems():
        # file path
        sort_bam = star_bam_path[:-3] + 'sort.bam'
        filter_bam = star_bam_path[:-3] + 'filter.bam'
        # command
        sort_cmd = f'samtools sort -o {sort_bam} --threads 2 {star_bam_path}'
        filter_cmd = f'samtools view -b -h -q {mapq_threshold} -o {filter_bam} {sort_bam}'
        cleaning_cmd = f'rm -f {star_bam_path} {sort_bam}'
        command = ' && '.join([sort_cmd, filter_cmd, cleaning_cmd])
        records.append([uid, index_name, filter_bam])
        command_list.append(command)

    with open(output_dir / 'star_bam_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'bam_path'])
    record_df.to_csv(output_dir / 'star_bam_qc.records.csv', index=None)
    return
