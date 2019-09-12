import glob
import logging
import pathlib
from collections import defaultdict

import pandas as pd

from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def merge_lane(output_dir: str):
    output_dir = pathlib.Path(output_dir)
    demultiplex_records_df = pd.read_csv(output_dir / 'demultiplex.records.csv')

    # make file sets dict, each kv pair will have a command and corresponding to one fastq file
    # key is (uid, index_name, read_type)
    file_set_dict = defaultdict(list)
    for uid, sub_df in demultiplex_records_df.groupby('uid'):
        for _, (_, _, r1_path, r2_path) in sub_df.iterrows():
            index_name = r1_path.split('-')[-2]  # based on demultiplex
            file_set_dict[(uid, index_name, 'R1')].append(r1_path)
            index_name = r2_path.split('-')[-2]  # based on demultiplex
            file_set_dict[(uid, index_name, 'R2')].append(r2_path)

    records = []
    command_list = []
    for (uid, index_name, read_type), file_path_list in file_set_dict.items():
        output_path = f'{output_dir}/{uid}-{index_name}-{read_type}.raw.fq.gz'
        path_str = ' '.join(file_path_list)
        merge_cmd = f'gzip -cd {path_str} | gzip -c > {output_path} && rm -f {path_str}'
        records.append([uid, index_name, read_type, output_path])
        command_list.append(merge_cmd)

    with open(output_dir / 'merge_lane.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records, columns=['uid', 'index_name', 'read_type', 'fastq_path'])
    record_df.to_csv(output_dir / 'merge_lane.records.csv', index=None)
    return record_df, command_list
