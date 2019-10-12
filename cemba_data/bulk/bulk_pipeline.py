import pathlib

import pandas as pd

from .allc_to_bigwig import generate_bigwig
from .extract_allc import extract_strand_merged_cg
from .merge_allc import _merge_cell, _merge_cluster


def bulk_pipeline(
        output_dir_path,
        group_table_path,
        chrom_size_path,
        binarize_single_cell=False,
        merge_cpu=10,
        ignore_names=None,
        max_cell_group=999,
        cg_context='CGN',
        bigwig_context=None,
        bigwig_binsize=50):
    output_dir = pathlib.Path(output_dir_path)
    output_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)

    group_table = pd.read_csv(group_table_path, index_col=0).fillna('NOT_A_CLUSTER')
    group_table_path = output_dir / 'GROUP_TABLE.csv'
    group_table.to_csv(group_table_path)

    # merge ALLC
    # step 1: merge cell ALLC
    command_file_list = []
    cell_group_allc_map, command_file_path = _merge_cell(group_table_path,
                                                         output_dir_path,
                                                         chrom_size_path,
                                                         binarize=binarize_single_cell,
                                                         cpu=merge_cpu,
                                                         max_cell_group=max_cell_group)
    command_file_list.append(command_file_path)

    # step 2: merge cluster ALLC
    if ignore_names is None:
        ignore_names = ['NOT_A_CLUSTER']
    else:
        ignore_names.append('NOT_A_CLUSTER')

    command_file_path = _merge_cluster(output_dir_path,
                                       cell_group_allc_map,
                                       ignore_names,
                                       chrom_size_path,
                                       cpu=merge_cpu)
    command_file_list.append(command_file_path)

    # step3: extract mCG ALLC
    command_file_path = extract_strand_merged_cg(output_dir_path,
                                                 chrom_size_path,
                                                 mc_context=cg_context,
                                                 cpu=1)
    command_file_list.append(command_file_path)

    # step4: bigwig
    if bigwig_context is not None:
        command_file_path = generate_bigwig(output_dir_path,
                                            chrom_size_path,
                                            cpu=1,
                                            mc_contexts_list=bigwig_context,
                                            bin_size=bigwig_binsize)
        command_file_list.append(command_file_path)
    return
