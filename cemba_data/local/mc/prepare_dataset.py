import json
import pathlib
import glob
import xarray as xr
import pandas as pd
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed


def generate_dataset(allc_files, out_dir, region_bed_path, region_name,
                     context_pattern, genome_size_path, max_cov_cutoff, dataset_name):
    """
    Batch version of map-to-region, a helper function that generating command.json for a whole dataset's ALLCs.
    The total number of region count bed files are len(region_name) * len(context_pattern)

    Parameters
    ----------
    allc_files
        List of ALLC files, accept wildcard
    out_dir
        Out_dir for all region count tables
    region_bed_path
        list of reference bed file paths
    region_name
        list of reference bed file names, corresponding to region_bed_path
    context_pattern
        list of mC context
    genome_size_path
        UCSC chrom.size like file, see instruction here: https://www.biostars.org/p/173963/
        Also can be downloaded from UCSC for most of the common genomes
    max_cov_cutoff
        maximum cov cutoff
    Returns
    -------

    """
    out_dir = pathlib.Path(out_dir).absolute()
    out_dir.mkdir(parents=True, exist_ok=True)

    # prepare out_dir
    cmd_list = []

    if '*' in allc_files:
        allc_files = glob.glob(allc_files)
    else:
        with open(allc_files) as f:
            allc_files = [line.strip() for line in f]

    region_paths = ' '.join(region_bed_path)
    region_names = ' '.join(region_name)
    for f in allc_files:
        f = pathlib.Path(f)
        out_path_prefix = str(out_dir / f.stem)
        _context_pattern = ' '.join(context_pattern)
        cmd = f'yap map-to-region --allc_path {f} ' \
              f'--out_path_prefix {out_path_prefix} ' \
              f'--region_bed_path {region_paths} ' \
              f'--region_name {region_names} ' \
              f'--genome_size_path {genome_size_path} ' \
              f'--context_pattern {_context_pattern} ' \
              f'--max_cov_cutoff {max_cov_cutoff}'
        command_dict = {
            'command': cmd,
            '-pe smp': 1,  # cpu for each command
        }
        cmd_list.append(command_dict)

    cmd_json_path = out_dir / 'map-to-region.command.json'
    with open(cmd_json_path, 'w') as f:
        json.dump(cmd_list, f)

    assemble_command = f'yap assemble-dataset --out_dir {out_dir}' \
                       f'--region_bed_path {region_paths}' \
                       f'--region_name {region_names}' \
                       f'--dataset_name {dataset_name}' \
                       f'--thread 5'
    assemble_json_path = out_dir / 'assemble.command.json'
    with open(assemble_json_path, 'w') as f:
        json.dump([{
            'command': assemble_command,
            '-pe smp': 10
        }], f)

    # submit master
    command_paths = ' '.join([cmd_json_path, assemble_json_path])
    qsub_command = f'yap qsub --working_dir {out_dir} ' \
                   f'--project_name generate-dataset ' \
                   f'--command_file_path {command_paths} ' \
                   f'--total_cpu 60 10'

    print(f"""
            The command file for generate-dataset is prepared
            ---------------------------------------------------------------------------------------
            - Output directory: {out_dir}
            - Yap command list: {command_paths}

            To run the command list using qsub, use "yap qsub" like this:

            {qsub_command}

            Modify the qsub parameters if you need. See "yap qsub -h" for help.
            """)

    return str(cmd_json_path)


def read_count_table(path):
    count_table = pd.read_table(path, index_col=3, header=None, na_values='.',
                                names=['chrom', 'start', 'end', 'id', 'mc', 'cov']).fillna(0)
    return count_table['mc'].astype(np.uint16).values, count_table['cov'].astype(np.uint16).values


def assemble_dataset(out_dir, region_bed_path, region_name, dataset_name,
                     thread=5):
    out_dir = pathlib.Path(out_dir).absolute()
    bed_paths = list(out_dir.glob('*count_table.bed.gz'))
    region_meta_dfs = {region_name: pd.read_table(region_bed, header=None, index_col=3,
                                                  names=['chrom', 'start', 'end', region_name])
                       for region_name, region_bed in zip(region_name, region_bed_path)}

    # Note: this part subject to the name pattern in generate_dataset
    records = []
    for path in bed_paths:
        id_part, _, type_part, *_ = (path.name.split('.'))
        cell = id_part.lstrip('allc_')
        region_type, mc_type = type_part.split('_')
        records.append([cell, region_type, mc_type, str(path)])
    count_table_df = pd.DataFrame(records,
                                  columns=['cell', 'region_type', 'mc_type', 'path']) \
        .set_index('cell')

    cell_index = count_table_df.index.unique()
    cell_count = cell_index.size
    mc_index = pd.Index(count_table_df.mc_type.unique())
    mc_type_count = mc_index.size

    da_dict = {}
    for region_type, _sub_df in count_table_df.groupby('region_type'):
        # initialize an empty ndarray with shape [cell, region, mc_type, count_type]
        region_meta = region_meta_dfs[region_type]
        region_count = region_meta.shape[0]
        region_data = np.zeros(shape=(cell_count, region_count, mc_type_count, 2),
                               dtype=np.uint16)

        for mc_i, mc_type in enumerate(mc_index):
            # make sure mc_type in order
            sub_df = _sub_df[_sub_df['mc_type'] == mc_type]
            # make sure cell in order
            sub_df = sub_df.reindex(cell_index)

            # multi-thread loading
            with ThreadPoolExecutor(max_workers=thread) as executor:
                future_data = {}
                for cell_i, path in enumerate(sub_df.path):
                    future_data[executor.submit(read_count_table, path)] = cell_i
                for future in as_completed(future_data):
                    cell_i = future_data[future]
                    mc, cov = future.result()
                    region_data[cell_i, :, mc_i, 0] = mc
                    region_data[cell_i, :, mc_i, 1] = cov
        region_dataarray = xr.DataArray(region_data,
                                        coords=[cell_index,
                                                region_meta.index,
                                                mc_index,
                                                pd.Index(['mc', 'cov'])],
                                        dims=['cell', region_type, 'mc_type', 'count_type'])
        da_dict[region_type + '_da'] = region_dataarray
    total_ds = xr.Dataset(da_dict)
    ds_out_path = out_dir / f'{dataset_name}.mcds'
    total_ds.to_netcdf(path=str(ds_out_path))
    return
