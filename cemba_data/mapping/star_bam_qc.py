import pathlib

import pandas as pd

from .utilities import get_configuration


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
        cleaning_cmd = ''  # f'rm -f {star_bam_path} {sort_bam}'
        command = ' && '.join([sort_cmd, filter_cmd, cleaning_cmd])
        records.append([uid, index_name, filter_bam])
        command_list.append(command)

    with open(output_dir / 'star_bam_qc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'bam_path'])
    record_df.to_csv(output_dir / 'star_bam_qc.records.csv', index=None)
    return record_df, command_list

