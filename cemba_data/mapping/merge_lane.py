import glob
import logging
import pathlib
import subprocess
from collections import defaultdict

import pandas as pd

from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def merge_lane(output_dir: str, config: str):
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)
    demultiplex_records_df = pd.read_csv(output_dir / 'demultiplex.records.csv')
    cutadapt_cores = config['fastqTrim']['cutadapt_cores']

    # make file sets dict, each kv pair will have a command and corresponding to one fastq file
    # key is (uid, index_name, read_type)
    file_set_dict = defaultdict(list)
    for uid, sub_df in demultiplex_records_df.groupby('uid'):
        for _, (_, _, r1_path_pattern, r2_path_pattern) in sub_df.iterrows():
            r1_path_list = glob.glob(r1_path_pattern)
            r2_path_list = glob.glob(r2_path_pattern)
            for path in r1_path_list:
                index_name = path.split('-')[-2]  # based on demultiplex
                file_set_dict[(uid, index_name, 'R1')].append(path)
            for path in r2_path_list:
                index_name = path.split('-')[-2]  # based on demultiplex
                file_set_dict[(uid, index_name, 'R2')].append(path)

    records = []
    command_list = []
    for (uid, index_name, read_type), file_path_list in file_set_dict.items():
        output_path = f'{output_dir}/{uid}-{index_name}-{read_type}.raw.fq.gz'
        path_str = ' '.join(file_path_list)
        merge_cmd = f'pigz -cd -p {cutadapt_cores} {path_str} | pigz -c > {output_path} && rm -f {path_str}'
        records.append([uid, index_name, read_type, output_path])
        command_list.append(merge_cmd)

    with open(output_dir / 'merge_lane.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records, columns=['uid', 'index_name', 'read_type', 'fastq_path'])
    record_df.to_csv(output_dir / 'merge_lane.records.csv', index=None)
    return


def merge_lane_runner(command):
    """
    merge fastq command wrapper, run command and parse results
    """
    try:
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       encoding='utf8',
                       shell=True,
                       check=True)
    except subprocess.CalledProcessError as e:
        log.error("Pipeline break, FASTQ merge lane ERROR! Command was:")
        log.error(command)
        log.error(e.stdout)
        log.error(e.stderr)
        raise e
    return
