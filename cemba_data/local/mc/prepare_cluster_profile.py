"""
Input:
- a dict json file. key: cluster ID, value: list of cell IDs
- metadata contain file path
- output dir

1. gather all file paths based on metadata
2. merge allc
3. extract CG strand merge allc
4. bigwig for CG and CH
5. file path summary file
"""
import pandas as pd
import pathlib
import json


def prepare_merge_allc(output_dir, group_table, chrom_size_file, skip_n_cell_cutoff=0,
                       cpu=10, h_vmem='5G', single_mode='copy'):
    """
    Prepare qsub command.json from a group table, following this format:
        1st column: the unique input file id
        2nd column: the group label
        3rd column: the absolute input file path
        Any other columns are ignored.
    This will generate all necessary file in the output_dir and can be directly used in qsub.

    Parameters
    ----------
    single_mode
        copy for copy the file, ln for create soft link
    output_dir
        Path to put output files, merged ALLC file will also save in this directory. Will create one if not exist.
    group_table
        Path or dataframe that contain group and input ALLC path information, first row must be column names
    chrom_size_file
        Path to UCSC chrom size file, chromosome names must agree with all the ALLC files.
    skip_n_cell_cutoff
        Skip small clusters with cell number <= cutoff, default is 0, which means not skip anything.
    cpu
        Number of CPUs to use in each merge command
    h_vmem
        MEM requirement for each merge command

    Returns
    -------

    """
    cpu = min(cpu, 20)

    if single_mode.lower() in ['copy', 'cp']:
        copy = True
    elif single_mode.lower() in ['ln', 'symlink']:
        copy = False
    else:
        raise ValueError(f"single_mode can only be in ['copy', 'cp', 'ln', 'symlink'], got {single_mode}")

    if isinstance(group_table, str):
        group_table = pd.read_csv(group_table, sep='\t', header=None)
    group_table.columns = ['uid', 'group', 'allc_path']
    group_table.set_index('uid', inplace=True)
    if group_table.index.duplicated().sum() != 0:
        raise ValueError('Input file ID in the first column of group_table is not unique.')

    path_dict = {}
    allc_list_dict = {}
    cmd_list = []
    for cluster, n_cell in group_table['group'].value_counts().iteritems():
        sub_df = group_table[group_table['group'] == cluster]
        n_cell = sub_df.shape[0]
        if n_cell <= skip_n_cell_cutoff:
            continue

        unique_file_path = output_dir / f'{cluster}.allc.tsv.gz'
        unique_cell_list_path = output_dir / f'{cluster}.allc_path.txt'

        path_dict[cluster] = str(unique_file_path)
        allc_list_dict[cluster] = sub_df['allc_path'].tolist()
        with open(unique_cell_list_path, 'w') as f:
            f.write('\n'.join(allc_list_dict[cluster]))

        if n_cell > 1:
            cmd = f'yap merge-allc --allc_paths {unique_cell_list_path} ' \
                  f'--out_path {unique_file_path} --chrom_size_file {chrom_size_file} --cpu {cpu}'
            cmd_dict = {'command': cmd, 'pe smp': int(cpu * 1.5), 'l h_vmem': h_vmem}
        else:
            # for cluster with only one cell (if possible), we just copy/symlink them to maintain structure
            if copy:
                cmd = f'cp {sub_df["allc_path"].iloc[0]} {unique_file_path};' \
                      f'cp {sub_df["allc_path"].iloc[0]}.tbi {unique_file_path}.tbi'
            else:
                cmd = f'ln -s {sub_df["allc_path"].iloc[0]} {unique_file_path};' \
                      f'cp {sub_df["allc_path"].iloc[0]}.tbi {unique_file_path}.tbi'
            cmd_dict = {'command': cmd, 'pe smp': 1, 'l h_vmem': '1G'}
        cmd_list.append(cmd_dict)

    with open(output_dir / 'GroupName_to_AllcPath.json', 'w') as f:
        json.dump(path_dict, f)
    with open(output_dir / 'command.json', 'w') as f:
        json.dump(cmd_list, f)
    return


def _generate_merge_strategy(cluster_table, min_group, keep_unique_cluster=True):
    """
    Generate merge strategy for multilayer cluster assignment. The cluster assignment do not need to be hierarchical.
    Through hierarchical form is better.
    1. cells that assigned same cluster among all columns are merged into cell groups, this is stage 0
    2. after stage 0, clusters are merged by cell groups in a minimum way, this is other stages

    """
    dropped_cells = []  # list of cell_id
    dropped_clusters = []  # (column, cluster_name)

    # make group dict, cell in each group have all-level same cluster assignment
    cell_group_dict = {}
    for group_id, (group_pattern, sub_df) in enumerate(cluster_table.groupby(cluster_table.columns.tolist())):
        if sub_df.shape[0] < min_group:
            dropped_cells += sub_df.index.tolist()
            # if a group have too less cell, this group is dropped
            continue
        cell_group_dict[group_id] = {
            'pattern': group_pattern,
            'cell_ids': sub_df.index,
            'cell_count': sub_df.shape[0],
            'out_path': None  # will be assign in _batch_merge_allc
        }

    # make cluster dict, cluster is merged from 1 or more groups
    cluster_dict = {}
    existing_group_pattern = {}
    real_cluster_name_map = {}
    for column_num, column in enumerate(cluster_table.columns):
        for cluster_name, sub_df in cluster_table.groupby(column):
            cluster_uid = f'{column}_{cluster_name}'

            # [real_column_dir, real_cluster_uid, original column, original cluster_name]
            # the 1st 2 items are used for make real path, last 2 items are used for sync original cluster table
            real_cluster_name_map[cluster_uid] = [column, cluster_uid, column, cluster_name]
            contain_groups = tuple(set([group_id for group_id, group_id_dict in cell_group_dict.items()
                                        if group_id_dict['pattern'][column_num] == cluster_name]))
            if len(contain_groups) == 0:
                dropped_clusters.append((column, cluster_name))
                # this cluster have no group, which happens to small cluster whose group dropped by min_group
                continue
            if keep_unique_cluster and (contain_groups in existing_group_pattern):
                # always keep cluster exists in later columns, delete the previous one.
                # assume that later columns are always in higher level.
                # this do not check the cell, but check on group level.
                del cluster_dict[existing_group_pattern[contain_groups]]

                # its ok to delete duplicated cluster here, but in the end,
                # we need to give out a file list to make it clear
                real_cluster_name_map[existing_group_pattern[contain_groups]][0] = column
                real_cluster_name_map[existing_group_pattern[contain_groups]][1] = cluster_uid

                # in case of chain deletion, search all k:v to replace v to the final cluster_uid
                for k, v in real_cluster_name_map.items():
                    if v[1] == existing_group_pattern[contain_groups]:
                        real_cluster_name_map[k][0] = column
                        real_cluster_name_map[k][1] = cluster_uid

            existing_group_pattern[contain_groups] = cluster_uid

            cluster_dict[cluster_uid] = {
                'groups': contain_groups,
                'cell_count': sub_df.shape[0],
                'out_path': None,  # will be assign in _batch_merge_allc
                'cluster_level': column
            }
    return cell_group_dict, cluster_dict, dropped_cells, dropped_clusters, real_cluster_name_map


def _batch_merge_allc(cluster_table, cell_path_series,
                      out_dir, min_group, cpu, chrom_size_file, bin_length):
    # raise NotImplementedError
    """
    Batch merge ALLC function, also accept multilayer cluster assignment.
    """
    # making directories
    out_dir = pathlib.Path(out_dir).absolute()
    out_dir.mkdir(exist_ok=False)
    group_dir = out_dir / 'merge_tmp_cell_group'
    group_dir.mkdir()
    for column in cluster_table.columns:
        column_dir = out_dir / column
        column_dir.mkdir()

    # generate merge dict
    cell_group_dict, cluster_dict, dropped_cells, dropped_clusters, real_cluster_name_map = \
        _generate_merge_strategy(cluster_table, min_group)

    # generate cell_group merge commands
    records = []
    for group_id, group_dict in cell_group_dict.items():
        cell_paths = cell_path_series.reindex(group_dict['cell_ids'])
        cell_id_list_path = group_dir / f'{group_id}.cell_list'
        with cell_id_list_path.open('w') as f:
            for path in cell_paths:
                f.write(str(path) + '\n')
        group_allc_out_path = group_dir / f'{group_id}.allc.tsv.gz'
        cell_group_dict[group_id]['out_path'] = group_allc_out_path
        cmd = f'yap merge-allc --allc_paths {cell_id_list_path} ' \
              f'--out_path {group_allc_out_path} --cpu {cpu} ' \
              f'--chrom_size_file {chrom_size_file} ' \
              f'--bin_length {bin_length}'
        memory_gbs = 5
        cmd_dict = {
            'command': cmd,
            'pe smp': int(cpu),
            'l h_vmem': f'{memory_gbs}G'
        }
        records.append(cmd_dict)
    group_command_path = out_dir / 'group_merge.command.json'
    with group_command_path.open('w') as f:
        json.dump(records, f)

    # generate cluster merge commands
    records = []
    for cluster_id, cluster_id_dict in cluster_dict.items():
        cluster_column = cluster_id_dict['cluster_level']
        cluster_column_dir = out_dir / cluster_column
        cluster_allc_out_path = cluster_column_dir / f'{cluster_id}.allc.tsv.gz'
        cluster_dict[cluster_id]['out_path'] = cluster_allc_out_path
        group_out_paths = [cell_group_dict[group_id]['out_path'] for group_id in cluster_id_dict['groups']]
        if len(group_out_paths) == 1:
            # no need to merge, copy the file instead,
            # through increase space usage, but most straightforward
            cmd1 = f'cp {group_out_paths[0]} {cluster_allc_out_path}'
            cmd2 = f'cp {group_out_paths[0]}.tbi {cluster_allc_out_path}.tbi'
            cmd_dicts = [{
                'command': cmd1,
                'pe smp': 10,  # TODO: figure out a better way to limit IO bound jobs
                'l h_vmem': '3G'
            }, {
                'command': cmd2,
                'pe smp': 10,  # TODO: figure out a better way to limit IO bound jobs
                'l h_vmem': '3G'
            }]
        else:
            group_id_list_path = cluster_column_dir / f'{cluster_id}.group_list'
            with group_id_list_path.open('w') as f:
                for path in group_out_paths:
                    f.write(str(path) + '\n')
            cmd = f'yap merge-allc --allc_paths {group_id_list_path} ' \
                  f'--out_path {cluster_allc_out_path} --cpu {cpu} ' \
                  f'--chrom_size_file {chrom_size_file} ' \
                  f'--bin_length {bin_length}'
            memory_gbs = 5
            cmd_dicts = [{
                'command': cmd,
                'pe smp': int(cpu * 1.1),
                'l h_vmem': f'{memory_gbs}G'
            }]
        records += cmd_dicts
    cluster_command_path = out_dir / 'cluster_merge.command.json'
    with cluster_command_path.open('w') as f:
        json.dump(records, f)

    # if keep_unique_cluster = True in _generate_merge_strategy,
    # some lower level duplicated cluster are dropped
    # here provide a final cluster path table follow original column and cluster name
    # where no cluster are dropped
    records = []
    for cluster_uid, (real_column, real_cluster_uid, column, cluster_name) in real_cluster_name_map.items():
        real_file_path = out_dir / real_column / f'{real_cluster_uid}.allc.tsv.gz'
        records.append([column, cluster_name, real_file_path])
    final_cluster_path_df = pd.DataFrame(records,
                                         columns=['column', 'cluster_name', 'cluster_allc_path'])
    final_cluster_path_df.to_csv(out_dir / 'final_cluster_path.tsv', index=None, sep='\t')
    return


def _batch_allc_profile(out_dir):
    records = []
    final_cluster_path_df = pd.read_table(out_dir / 'final_cluster_path.tsv', header=0)
    allc_files = final_cluster_path_df['cluster_allc_path'].drop_duplicates()
    for allc_file in allc_files:
        cmd = f'yap allc-profile --allc_path {allc_file} --out_path {allc_file}.profile'
        cmd_dict = {
            'command': cmd,
            'pe smp': 1,
            'l h_vmem': '3G'
        }
        records.append(cmd_dict)
    command_path = out_dir / 'allc_profile.command.json'
    with command_path.open('w') as f:
        json.dump(records, f)
    return


def _batch_allc_to_bigwig(out_dir, chrom_size_path, mc_contexts):
    records = []
    final_cluster_path_df = pd.read_table(out_dir / 'final_cluster_path.tsv', header=0)
    allc_files = final_cluster_path_df['cluster_allc_path'].drop_duplicates()
    for mc_context in mc_contexts:
        for allc_file in allc_files:
            out_path = str(allc_file).rstrip('tsv.gz') + f'{mc_context}.bw'
            cmd = f'yap allc-to-bigwig --allc_path {allc_file} ' \
                  f'--out_path {out_path} --chrom_size {chrom_size_path} --mc_type {mc_context}'
            cmd_dict = {
                'command': cmd,
                'pe smp': 1,
                'l h_vmem': '3G'
            }
            records.append(cmd_dict)
    command_path = out_dir / 'allc_bigwig.command.json'
    with command_path.open('w') as f:
        json.dump(records, f)
    return


def _batch_extract_mc(out_dir, mc_contexts, merge_strand):
    if isinstance(mc_contexts, str):
        mc_contexts = [mc_contexts]
    records = []
    final_cluster_path_df = pd.read_table(out_dir / 'final_cluster_path.tsv', header=0)
    allc_files = final_cluster_path_df['cluster_allc_path'].drop_duplicates()
    for mc_context in mc_contexts:
        for allc_file in allc_files:
            out_path = str(allc_file).rstrip('tsv.gz') + f'{mc_context}.tsv.gz'
            cmd = f'yap allc-extract --allc_path {allc_file} ' \
                  f'--out_path {out_path} --merge_strand {merge_strand} --mc_context {mc_context}'
            cmd_dict = {
                'command': cmd,
                'pe smp': 1,
                'l h_vmem': '3G'
            }
            records.append(cmd_dict)
    command_path = out_dir / 'allc_extract.command.json'
    with command_path.open('w') as f:
        json.dump(records, f)
    return


def cluster_merge_pipeline(cluster_table_path, allc_path_file, out_dir,
                           chrom_size_path, bin_length=1000000,
                           bigwig_contexts=('CGN', 'CHN'),
                           extract_contexts=('CGN',), merge_strand=True,
                           min_group=10, merge_allc_cpu=20, total_cpu=80):
    """
    Merge allc and generate bigwig pipeline

    Parameters
    ----------
    cluster_table_path
        Path to the cluster table WITH 1 HEADER ROW. The first column must be cell id. The other columns
        are different levels of cluster assignments. From left to right, sub to major.
    allc_path_file
        Path to the ALLC path table WITHOUT HEADER. The first column must be cell id. The second column
        is the ALLC file path for each cell.
    out_dir
        Output directory, must not exist.
    chrom_size_path
        Path to UCSC chrom size file, used for guide merge and bigwig generation.
    bin_length
        Length of the chrom bin size when do parallel merging. The larger the more MEM usage.
        The default value is usually fine.
    bigwig_contexts
        mC contexts used for generate bigwig
    extract_contexts
        mC contexts used for extract sub ALLC
    merge_strand
        Whether to merge adjacent CG site methylation info
    min_group
        Minimum number of cells for a cell group to be considered.
    merge_allc_cpu
        CPU used in each merge allc job.
    total_cpu
        Total CPU used in qsub queue
    Returns
    -------

    """

    cluster_table = pd.read_table(cluster_table_path, index_col=0, header=0)
    if cluster_table.columns.duplicated().sum() != 0:
        raise ValueError('Cluster table have duplicated columns')
    if cluster_table.columns.duplicated().sum() != 0:
        raise ValueError('Cluster table have duplicated cell_ids')

    cell_path_series = pd.read_table(allc_path_file, index_col=0, header=None, squeeze=True)

    out_dir = pathlib.Path(out_dir).absolute()
    _batch_merge_allc(cluster_table, cell_path_series=cell_path_series,
                      out_dir=out_dir, min_group=min_group, cpu=merge_allc_cpu,
                      chrom_size_file=chrom_size_path, bin_length=bin_length)
    _batch_allc_profile(out_dir=out_dir)
    _batch_allc_to_bigwig(out_dir, chrom_size_path, mc_contexts=bigwig_contexts)
    _batch_extract_mc(out_dir, mc_contexts=extract_contexts, merge_strand=merge_strand)

    command_path_list = [
        out_dir / 'group_merge.command.json',
        out_dir / 'cluster_merge.command.json',
        out_dir / 'allc_profile.command.json',
        out_dir / 'allc_bigwig.command.json',
        out_dir / 'allc_extract.command.json'
    ]

    # submit master of master
    command_path = '"' + '" "'.join(map(str, command_path_list)) + '"'
    qsub_command = f'yap qsub --working_dir {out_dir} ' \
                   f'--project_name merge ' \
                   f'--command_file_path {command_path} ' \
                   f'--total_cpu {total_cpu} ' \
                   f'--total_mem 1000 '

    print(f"""
        The command file for prepare-cluster-profile is prepared
        ---------------------------------------------------------------------------------------
        - Output directory: {out_dir}
        - Yap command list: {command_path}

        To run the command list using qsub, use "yap qsub" like this:

        {qsub_command}

        Modify the qsub parameters if you need. See "yap qsub -h" for help.
        """)

    return
