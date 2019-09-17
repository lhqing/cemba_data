import json
import pathlib
from collections import defaultdict

import pandas as pd


def prepare_merge_first_level(cell_cluster_table,
                              output_dir,
                              cluster_col_name,
                              chrom_size_path,
                              min_cell=0,
                              max_cell=1000,
                              seed=1,
                              cpu=20):
    """The first column is merged from single cell ALLC files"""
    qsub_dir = output_dir / 'qsub/merge_allc'
    qsub_dir.mkdir(exist_ok=True)
    col_output_dir = output_dir / cluster_col_name
    col_output_dir.mkdir(exist_ok=True)

    command_records = []
    output_records = {}
    for cluster, sub_df in cell_cluster_table.groupby(cluster_col_name):
        allc_series = sub_df.iloc[:, 0]  # the first column is always ALLC path
        if allc_series.size < min_cell:  # skip very small clusters
            continue
        if allc_series.size > max_cell:
            allc_series = allc_series.sample(max_cell, random_state=seed)
        output_path = col_output_dir / f'{cluster}.allc.tsv.gz'
        output_records[cluster] = str(output_path)

        cell_list_path = col_output_dir / f'{cluster}.allc_paths.txt'
        with open(cell_list_path, 'w') as f:
            for path in allc_series.values:
                f.write(str(path) + '\n')
        cmd = f'allcools merge-allc ' \
              f'--allc_paths {cell_list_path} ' \
              f'--output_path {output_path} ' \
              f'--chrom_size_path {chrom_size_path} ' \
              f'--cpu {cpu}'
        command_records.append(cmd)
    with open(qsub_dir / f'merge_{cluster_col_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(col_output_dir / f'merge_{cluster_col_name}.records.json', 'w') as f:
        json.dump(output_records, f)
    return output_records


def merge_other_levels(output_dir,
                       cell_cluster_table,
                       cur_col_name,
                       source_col_name,
                       source_allc_dict,
                       chrom_size_path,
                       cpu):
    """The other columns are merged from first column ALLC files"""
    qsub_dir = output_dir / 'qsub/merge_allc'
    col_output_dir = output_dir / cur_col_name
    col_output_dir.mkdir(exist_ok=True)

    cluster_map = defaultdict(list)

    cluster_df = cell_cluster_table[[source_col_name, cur_col_name]]
    for _, (source_col_cluster, cur_cluster) in cluster_df.drop_duplicates().iterrows():
        try:
            source_path = source_allc_dict[source_col_cluster]
        except KeyError:
            # KeyError means the source_col_cluster is skipped in last level
            # if the only source_col_cluster is skipped, then this cluster will be skipped too.
            continue
        cluster_map[cur_cluster].append(source_path)

    command_records = []
    output_records = {}
    for cluster, allc_list in cluster_map.items():
        output_path = f'{col_output_dir}/{cluster}.allc.tsv.gz'
        output_records[cluster] = str(output_path)

        if len(allc_list) == 1:
            cmd = f'ln -s {allc_list[0]} {output_path} && ln -s {allc_list[0]}.tbi {output_path}.tbi'
            # only put real cluster
        else:
            allc_list_str = ' '.join(allc_list)
            cmd = f'allcools merge-allc ' \
                  f'--allc_paths {allc_list_str} ' \
                  f'--output_path {output_path} ' \
                  f'--chrom_size_path {chrom_size_path} ' \
                  f'--cpu {cpu}'
        allc_list_path = col_output_dir / f'{cluster}.allc_paths.txt'
        with open(allc_list_path, 'w') as f:
            for path in allc_list:
                f.write(str(path) + '\n')
        command_records.append(cmd)

    with open(qsub_dir / f'merge_{cur_col_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(col_output_dir / f'merge_{cur_col_name}.records.json', 'w') as f:
        json.dump(output_records, f)
    return output_records


def extract_cg(output_dir, sub_dir_name, chrom_size_path, cpu, mc_contexts='CGN'):
    qsub_dir = output_dir / 'qsub/extract_cg'
    qsub_dir.mkdir(exist_ok=True)

    _output_dir = output_dir / sub_dir_name
    output_records_path = _output_dir / f'merge_{sub_dir_name}.records.json'
    with open(output_records_path) as f:
        allc_dict = json.load(f)
    command_records = []
    output_records = {}
    for cluster, allc_path in allc_dict.items():
        output_prefix = _output_dir / cluster
        output_path = _output_dir / (cluster + f'.{mc_contexts}-Merge.allc.tsv.gz')
        output_records[cluster] = output_path
        if allc_path.is_symlink():
            # for these smaller files, use real command and actual computation is clearer
            allc_path = allc_path.resolve()
        command = f'allcools extract-allc ' \
                  f'--allc_path {allc_path} ' \
                  f'--output_prefix {output_prefix} ' \
                  f'--mc_contexts {mc_contexts} ' \
                  f'--chrom_size_path {chrom_size_path} ' \
                  f'--strandness merge ' \
                  f'--output_format allc ' \
                  f'--cpu {cpu}'
        command_records.append(command)
    with open(qsub_dir / f'extract_cg_{sub_dir_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(_output_dir / f'extract_cg_{sub_dir_name}.records.json', 'w') as f:
        json.dump(output_records, f)


def generate_bigwig(output_dir,
                    sub_dir_name,
                    chrom_size_path,
                    cpu,
                    mc_contexts_list,
                    bin_size):
    qsub_dir = output_dir / 'qsub/bigwig'
    qsub_dir.mkdir(exist_ok=True)

    _output_dir = output_dir / sub_dir_name
    output_records_path = _output_dir / f'merge_{sub_dir_name}.records.json'
    with open(output_records_path) as f:
        allc_dict = json.load(f)
    command_records = []
    output_records = {}
    for cluster, allc_path in allc_dict.items():
        output_prefix = _output_dir / cluster

        if allc_path.is_symlink():
            real_allc = allc_path.resolve()
            bw_files = real_allc.parent.glob(f'{cluster}*bw')
            for bw_file_path in bw_files:
                bw_suffix = '.'.join(bw_file_path.name.split('.')[-3:])
                command = f'ln -s {bw_file_path} {output_prefix}.{bw_suffix}'
                command_records.append(command)
        else:
            mc_contexts_str = ' '.join(mc_contexts_list)
            command = f'allcools allc-to-bigwig ' \
                      f'--allc_path {allc_path} ' \
                      f'--output_prefix {output_prefix} ' \
                      f'--chrom_size_path {chrom_size_path} ' \
                      f'--mc_contexts {mc_contexts_str} ' \
                      f'--bin_size {bin_size} ' \
                      f'--remove_additional_chrom ' \
                      f'--cpu {cpu}'
            command_records.append(command)
    with open(qsub_dir / f'bigwig_{sub_dir_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(_output_dir / f'bigwig_{sub_dir_name}.records.json', 'w') as f:
        json.dump(output_records, f)


def bulk_pipeline(
        output_path,
        cell_cluster_table_path,
        chrom_size_path,
        bw_context_list,
        bw_bin_size=50,
        max_cell=500,
        min_cell=10,
        seed=1,
        cpu=20,
        cg_context='CGN'):
    output_dir = pathlib.Path(output_path)
    output_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)

    cell_cluster_table = pd.read_csv(cell_cluster_table_path,
                                     sep='\t', index_col=0)
    first_col, *other_cols = cell_cluster_table.columns[1:]

    first_col_allc_dict = prepare_merge_first_level(
        cell_cluster_table=cell_cluster_table,
        output_dir=output_dir,
        cluster_col_name=first_col,
        chrom_size_path=chrom_size_path,
        min_cell=min_cell,
        max_cell=max_cell,
        seed=seed,
        cpu=cpu)

    # TODO submit to qsub

    for col_name in other_cols:
        merge_other_levels(
            output_dir=output_dir,
            cell_cluster_table=cell_cluster_table,
            source_col_name=first_col,
            cur_col_name=col_name,
            source_allc_dict=first_col_allc_dict,
            chrom_size_path=chrom_size_path,
            cpu=cpu)

    # TODO submit to qsub

    for col_name in cell_cluster_table.columns[1:]:
        extract_cg(output_dir=output_dir,
                   sub_dir_name=col_name,
                   chrom_size_path=chrom_size_path,
                   cpu=min(10, cpu),
                   mc_contexts=cg_context)

        generate_bigwig(output_dir=output_dir,
                        sub_dir_name=col_name,
                        chrom_size_path=chrom_size_path,
                        cpu=min(10, cpu),
                        mc_contexts_list=bw_context_list,
                        bin_size=bw_bin_size)

    # TODO submit to qsub
