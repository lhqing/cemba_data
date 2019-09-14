import logging
import pathlib

import pandas as pd

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def merge_bam(output_dir, record_name):
    """
    """
    output_dir = pathlib.Path(output_dir)
    bismark_records = pd.read_csv(output_dir / record_name,
                                  squeeze=True)

    # process bam
    records = []
    command_list = []
    for (uid, index_name), sub_df in bismark_records.groupby(['uid', 'index_name']):
        # file path
        r1_r2_bam_path_str = ' '.join(sub_df['bam_path'].tolist())
        output_bam = output_dir / f'{uid}_{index_name}.final.bam'
        # command
        command = f'samtools merge -f {output_bam} {r1_r2_bam_path_str} ' \
                  f'&& rm -f {r1_r2_bam_path_str}'
        records.append([uid, index_name, output_bam])
        command_list.append(command)

    with open(output_dir / 'final_bam.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'bam_path'])
    record_df.to_csv(output_dir / 'final_bam.records.csv', index=None)
    return record_df, command_list
