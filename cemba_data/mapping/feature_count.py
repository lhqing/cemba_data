import pathlib
import pandas as pd
import numpy as np
import subprocess
import shlex
from concurrent.futures import ProcessPoolExecutor, as_completed

from cemba_data.mapping.utilities import get_configuration

def _execute_one_feature_count(record):
    try:
        subprocess.run(shlex.split(record['cmd']), check=True)
    except Exception as e:
        print('This featureCount command generate error:')
        print(record['cmd'])
        raise e
    return


def _read_one_feature_count_output(record):
    out_path = record['out_path']
    summary_path = out_path + '.summary'

    summary_df = pd.read_table(summary_path, index_col=0,
                               names=record['samples'], skiprows=1)

    count_df = pd.read_table(out_path, comment='#', index_col=0)
    count_df.columns = count_df.columns[:5].tolist() + record['samples']
    count_df = count_df.iloc[:, 5:].astype(np.uint32)
    count_df.index.name = 'gene'
    return count_df, summary_df


def _assemble_feature_count_output(records):
    count_dfs = []
    summary_dfs = []
    for record in records:
        count_df, summary_df = _read_one_feature_count_output(record)
        count_dfs.append(count_df)
        summary_dfs.append(summary_df)
    total_count = pd.concat(count_dfs, axis=1, sort=True)
    total_summary = pd.concat(summary_dfs, axis=1, sort=True).T
    return total_summary, total_count


def batch_feature_count(bam_table, out_prefix, gtf_path,
                        count_type='gene', id_type='gene_id',
                        cpu=2, chunksize=50):
    """
    Count RNA read using featureCount, return a pandas msgpack file

    Parameters
    ----------
    bam_table
        tab separated, first column is cell id, second column is file path, no header.
    out_prefix
    gtf_path
    count_type
    id_type
    cpu
    chunksize

    Returns
    -------

    """
    bam_dict = pd.read_csv(bam_table,
                           sep='\t',
                           index_col=0,
                           header=None,
                           squeeze=True).to_dict()

    parent_dir = pathlib.Path(out_prefix).parent
    if not parent_dir.exists():
        raise FileNotFoundError('Parent directory of out_prefix do not exist.')

    count = 0
    records = []
    sample_ids = []
    bam_paths = []

    for sample_id, bam_path in bam_dict.items():
        sample_ids.append(sample_id)
        bam_paths.append(bam_path)
        count += 1
        if count % chunksize == 0:
            chunk_id = count // chunksize
            bam_paths_str = ' '.join(bam_paths)
            out_path = out_prefix + f'.tmp.{count_type}.{id_type}.{chunk_id}.tsv'
            cmd = f'featureCounts -t {count_type} ' \
                  f'-g {id_type} -a {gtf_path} -o {str(out_path)} {bam_paths_str}'
            cmd_dict = {
                'cmd': cmd,
                'samples': sample_ids,
                'out_path': str(out_path)
            }
            records.append(cmd_dict)
            sample_ids = []
            bam_paths = []
    # last chunk:
    chunk_id = count // chunksize + 1  # prevent dup
    bam_paths_str = ' '.join(bam_paths)
    out_path = out_prefix + f'.tmp.{count_type}.{id_type}.{chunk_id}.tsv'
    cmd = f'featureCounts -t {count_type} ' \
          f'-g {id_type} -a {gtf_path} -o {str(out_path)} {bam_paths_str}'
    cmd_dict = {
        'cmd': cmd,
        'samples': sample_ids,
        'out_path': str(out_path)
    }
    records.append(cmd_dict)

    with ProcessPoolExecutor(cpu) as executor:
        future_dict = {executor.submit(_execute_one_feature_count,
                                       record=record): record_id
                       for record_id, record in enumerate(records)}

        for future in as_completed(future_dict):
            record_id = future_dict[future]
            future.result()  # try to get return status or the error will raise

    total_summary, total_count = _assemble_feature_count_output(records)

    output_store_path = out_prefix + '.h5'
    with pd.HDFStore(output_store_path, 'a') as store:
        store[f'{count_type}-{id_type}'] = total_count
        store[f'{count_type}-{id_type}-summary'] = total_summary

    # rm tmp files
    subprocess.run(f'rm -f {out_prefix}*tmp*', shell=True)
    return


def prepare_feature_count(output_dir, config):
    config = get_configuration(config)

    final_rna_bam_records = pd.read_csv(output_dir / 'select_rna_reads.records.csv')
    final_rna_bam_records.index = final_rna_bam_records['uid'] + '-' + final_rna_bam_records['index_name']
    bam_table_path = output_dir / 'bam_table_for_feature_count.tsv'
    output_feature_count_prefix = output_dir / 'featureCount'
    gtf_path = config['featureCount']['gtf_path']
    count_type = config['featureCount']['count_type']
    id_type = config['featureCount']['id_type']

    with open(bam_table_path, 'w') as f:
        for cell_id, bam_path in final_rna_bam_records['bam_path'].iteritems():
            f.write(f'{cell_id}\t{bam_path}\n')
    feature_count_command = f'yap-internal featurecount ' \
                            f'--bam_table {bam_table_path} ' \
                            f'--out_prefix {output_feature_count_prefix} ' \
                            f'--gtf_path {gtf_path} ' \
                            f'--count_type {count_type} ' \
                            f'--id_type {id_type} ' \
                            f'--cpu 10 ' \
                            f'--chunksize 50'
    with open(output_dir / 'feature_count.records.txt', 'w') as f:
        f.write(feature_count_command)
    with open(output_dir / 'feature_count.command.txt', 'w') as f:
        output_path = f'{output_feature_count_prefix}.h5'
        f.write(output_path)
    return output_path, feature_count_command

