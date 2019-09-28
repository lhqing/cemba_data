import pathlib
from ..mapping.utilities import command_runner
from .merge_allc import _merge_cell, _merge_cluster
from ..qsub import qsub


def bulk_pipeline(
        output_dir_path,
        group_table_path,
        chrom_size_path,
        binarize_single_cell=True,
        merge_cpu=10,
        total_cpu=100,
        ignore_names=None,
        mode='command_only',
        h_vmem='5G',
        max_cell_group=300):
    output_dir = pathlib.Path(output_dir_path)
    output_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)

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
    command_file_path = _merge_cluster(output_dir_path,
                                       cell_group_allc_map,
                                       ignore_names,
                                       chrom_size_path,
                                       cpu=merge_cpu)
    command_file_list.append(command_file_path)

    # merge ALLC runner
    for command_file_path in command_file_list:
        if mode == 'qsub':
            qsub(command_file_path=str(command_file_path),
                 working_dir=qsub_dir,
                 project_name='merge',
                 wait_until=None,
                 total_cpu=total_cpu,
                 total_mem=500,
                 force_redo=False,
                 qsub_global_parms=f'-pe smp {merge_cpu};-l h_vmem={h_vmem}')
        elif mode == 'command_only':
            pass
        elif mode == 'local':
            with open(command_file_path) as f:
                commands = [l.strip() for l in f.readlines()]
            total_cpu = min(30, total_cpu)
            command_runner(commands, cpu=total_cpu)
        else:
            raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')
    return
