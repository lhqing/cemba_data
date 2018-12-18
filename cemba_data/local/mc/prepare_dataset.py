import json
import pathlib
import glob
import xarray as xr
import pandas as pd
import numpy as np
import multiprocessing
import argparse


def batch_map_to_region(allc_files, out_dir, region_bed_path, region_name,
                        context_pattern, genome_size_path, max_cov_cutoff):
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

    for f in allc_files:
        if 'stats' in f.parent.name:
            continue
        out_path_prefix = str(out_dir / f.stem)
        region_paths = ' '.join(region_bed_path)
        region_names = ' '.join(region_name)
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

    cmd_json_path = out_dir / 'map-to-region_command.json'
    with open(cmd_json_path, 'w') as f:
        json.dump(cmd_list, f)
    qsub_command = f'yap qsub --working_dir {out_dir} ' \
                   f'--project_name map-to-region ' \
                   f'--command_file_path {cmd_json_path} ' \
                   f'--total_cpu 160 --submission_gap 2 --qstat_gap 30'

    print(f"""
    The command file for map-to-region is prepared
    ---------------------------------------------------------------------------------------
    - Output directory: {out_dir}
    - Yap mapping command list: {cmd_json_path}

    To run the command list using qsub, use "yap qsub" like this:

    {qsub_command}

    Modify the qsub parameters if you need. See "yap qsub -h" for help.
    """)

    return str(cmd_json_path)


def _read_count_table(file_path, region_name, cell_id):
    region_df = pd.read_table(file_path, na_values='.', index_col=region_name,
                              names=['chrom', 'start', 'end', region_name, 'mc', 'cov']).fillna(0)
    region_mc = region_df['mc'].astype(np.uint16)
    region_mc.rename(cell_id, inplace=True)
    region_cov = region_df['cov'].astype(np.uint16)
    region_cov.rename(cell_id, inplace=True)
    return region_mc, region_cov


def assemble_dataset(out_dir, dataset_name, cpu):
    records = []
    for f in out_dir.glob('**/*.count_table.bed.gz'):
        id_part, type_part, _, _, _ = f.name.split('.')
        _, uid, index_name = id_part.split('_')
        region_name, mc_type = type_part.split('_')
        records.append({'uid': uid,
                        'index_name': index_name,
                        'region_name': region_name,
                        'mc_type': mc_type,
                        'file_path': f})
    count_path_df = pd.DataFrame(records)

    region_da_dict = {}
    for region_name, region_df in count_path_df.groupby('region_name'):
        mc_das = []
        for mc_type, mc_type_df in region_df.groupby('mc_type'):
            pool = multiprocessing.Pool(cpu)
            # read bed files
            results = []
            for _, row in mc_type_df.iterrows():
                cell_id = f"{row['uid']}-{row['index_name']}"
                result = pool.apply_async(_read_count_table,
                                          kwds=dict(file_path=row['file_path'],
                                                    region_name=row['region_name'],
                                                    cell_id=cell_id))
                results.append(result)
            pool.close()
            pool.join()

            # make region data array
            region_das = []
            for result in results:
                region_mc, region_cov = result.get()
                region_mc_da = region_mc.to_xarray().expand_dims(['count_type', 'cell'])
                region_mc_da.coords['count_type'] = ['mc']
                region_mc_da.coords['cell'] = [region_cov.name]
                region_cov_da = region_cov.to_xarray().expand_dims(['count_type', 'cell'])
                region_cov_da.coords['count_type'] = ['cov']
                region_cov_da.coords['cell'] = [region_cov.name]
                region_das.append(xr.concat([region_mc_da, region_cov_da],
                                            dim='count_type'))
            region_total_da = xr.concat(region_das, dim='cell').expand_dims(['mc_type'])
            region_total_da.coords['mc_type'] = [mc_type]
            mc_das.append(region_total_da)
        region_da_dict[region_name + '_da'] = xr.concat(mc_das, dim='mc_type')
    total_dataset = xr.Dataset(region_da_dict)
    total_dataset.to_netcdf(path=f'{out_dir}/{dataset_name}.mcds')
    return

