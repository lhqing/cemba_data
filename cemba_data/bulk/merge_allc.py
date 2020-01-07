import json
import pathlib
from collections import defaultdict

import pandas as pd


def _merge_cell(group_table_path, output_dir_path, chrom_size_path, binarize, cpu, max_cell_group):
    output_dir = pathlib.Path(output_dir_path).absolute()
    qsub_dir = output_dir / 'qsub'

    group_table = pd.read_csv(group_table_path, index_col=0)

    if binarize:
        binarize_param = '--binarize'
    else:
        binarize_param = ''

    # cell level merge
    cell_merge_output_dir = output_dir / 'CELL_GROUP'
    cell_merge_output_dir.mkdir(exist_ok=True)

    command_records = {}
    col_names = group_table.columns[1:].tolist()
    file_map = {col_name: defaultdict(list) for col_name in col_names}
    for i, (cluster_combination, sub_df) in enumerate(group_table.groupby(col_names)):
        cell_number = sub_df.shape[0]
        if (max_cell_group is not None) and (cell_number > max_cell_group):
            sub_df = sub_df.sample(max_cell_group, random_state=1)
            cell_number = max_cell_group

        output_path = cell_merge_output_dir / f'{i}-{cell_number}.allc.tsv.gz'
        file_list_path = cell_merge_output_dir / f'{i}-{cell_number}.cell_list.txt'

        if isinstance(cluster_combination, str):
            cluster_combination = [cluster_combination]

        for col_name, cluster_name in zip(col_names, cluster_combination):
            file_map[col_name][cluster_name].append(str(output_path))

        # dump file list
        with open(file_list_path, 'w') as f:
            f.write('\n'.join(sub_df['AllcPath']))

        # command
        if cell_number == 1:
            input_path = sub_df['AllcPath'][0]
            cmd = f'cp {input_path} {output_path}; ' \
                  f'cp {input_path}.tbi {output_path}.tbi'
        else:
            input_path = file_list_path
            cmd = f'allcools merge-allc ' \
                  f'--allc_paths {input_path} ' \
                  f'--output_path {output_path} ' \
                  f'--chrom_size_path {chrom_size_path} ' \
                  f'--cpu {cpu} {binarize_param}'
        command_records[cmd] = cell_number

    command_path = qsub_dir / '_merge_cell_commands.txt'
    with open(command_path, 'w') as f:
        # write down commands order by cell number, larger number first
        for cmd in pd.Series(command_records).sort_values(ascending=False).index:
            f.write(cmd + '\n')

    with open(cell_merge_output_dir / '_file_map.json', 'w') as f:
        json.dump(file_map, f)
    return file_map, command_path


def _merge_cluster(output_dir_path, cell_group_allc_map, ignore_names, chrom_size_path, cpu):
    output_dir = pathlib.Path(output_dir_path).absolute()
    qsub_dir = output_dir / 'qsub'

    if ignore_names is None:
        ignore_names = []

    command_records = []
    for col, cluster_path_map in cell_group_allc_map.items():
        col_output_dir = output_dir / col
        col_output_dir.mkdir(exist_ok=True)
        file_path_records = []
        for cluster, cell_group_allc_list in cluster_path_map.items():
            if cluster in ignore_names:
                continue

            output_path = col_output_dir / f'{col}.{cluster}.allc.tsv.gz'

            total_cell_number = 0
            for cell_group_allc_path in cell_group_allc_list:
                # parse cell number from cell_group_allc file name, fix by _merge_cell_allc
                total_cell_number += int(pathlib.Path(cell_group_allc_path).name.split('.')[0].split('-')[-1])

            if len(cell_group_allc_list) == 1:
                cmd = f'cp {cell_group_allc_list[0]} {output_path}; ' \
                      f'cp {cell_group_allc_list[0]}.tbi {output_path}.tbi'
            else:
                input_path = col_output_dir / f'{col}.{cluster}.list.txt'
                with open(input_path, 'w') as f:
                    f.write('\n'.join(cell_group_allc_list))
                cmd = f'allcools merge-allc ' \
                      f'--allc_paths {input_path} ' \
                      f'--output_path {output_path} ' \
                      f'--chrom_size_path {chrom_size_path} ' \
                      f'--cpu {cpu}'
            command_records.append(cmd)
            file_path_records.append([col, cluster, output_path, total_cell_number])
        file_df = pd.DataFrame(file_path_records,
                               columns=['LevelName', 'ClusterName',
                                        'AllcPath', 'TotalCellCombined'])
        file_df.to_csv(col_output_dir / 'ALLC_table.csv', index=None)

    command_path = qsub_dir / '_merge_cluster_commands.txt'
    with open(command_path, 'w') as f:
        f.write('\n'.join(command_records))
    return command_path
