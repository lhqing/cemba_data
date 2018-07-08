import configparser
import os
import datetime
import json
import numpy as np
import pandas as pd
import h5py
import inspect
from numpy import uint32
from scipy.sparse import lil_matrix
from .qsub import Qsubmitter, qsub_config
from subprocess import run

ref_path_config = configparser.ConfigParser()
ref_path_config.read(os.path.dirname(__file__) + '/config_ref_path.ini')


def cur_time():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def batch_map_to_region(cell_ids, allc_files, out_dir, genome, dataset_name=None,
                        region_bed_path=None, region_name=None,
                        context_pattern=None, max_cov_cutoff=None,
                        remove_tmp=True, tmp_compression=False,
                        total_cpu=50, submission_gap=5,
                        qstat_gap=100, submit=False):
    # update qsub_config
    qsub_config['RUNNING_DEFAULT']['TOTAL_CPU'] = str(total_cpu)
    qsub_config['RUNNING_DEFAULT']['SUBMISSION_GAP'] = str(submission_gap)
    qsub_config['RUNNING_DEFAULT']['QSTST_GAP'] = str(qstat_gap)
    job_dir = qsub_config['QSUB_DEFAULT']['JOB_DIR']
    if dataset_name is None:  # if dataset name not provided, use time
        dataset_name = cur_time()
    project_name = f'map_to_region_{dataset_name}'
    working_job_dir = job_dir + f'/{project_name}'
    command_file_path = working_job_dir + '/command_list.json'
    try:
        os.mkdir(working_job_dir)
    except OSError:
        print('Note: The job has been submitted before.')
        print(f'If want to resubmit everything, delete or rename this directory: {working_job_dir}')
        pass

    # make refs
    genome = genome.upper()
    genome_size_path = ref_path_config[genome]['CHROM_SIZE']
    cell_level_region_set = []
    cell_level_region_name = []
    for _region_set, _region_name in zip(ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP'].split(' '),
                                         ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP_NAME'].split(' ')):
        cell_level_region_set.append(ref_path_config[genome][_region_set])
        cell_level_region_name.append(_region_name)
    if region_bed_path is None:
        region_bed_path = cell_level_region_set
    region_bed_path = ' '.join(region_bed_path)
    if region_name is None:
        region_name = cell_level_region_name
    region_name = ' '.join(region_name)
    if context_pattern is None:
        context_pattern = ref_path_config['DATA_SELECT']['ALLC_CONTEXT'].split(' ')
    context_pattern = ' '.join(context_pattern)

    # make command json
    command_dict_list = []
    if len(cell_ids) != len(allc_files):
        raise ValueError('Cell id and allc files have different length.')
    for cell, allc in zip(cell_ids, allc_files):
        command = f'yap map-to-region --allc_path {allc} ' \
                  f'--out_path_prefix {out_dir}/{cell} ' \
                  f'--region_bed_path {region_bed_path} ' \
                  f'--region_name {region_name} ' \
                  f'--genome_size_path {genome_size_path} ' \
                  f'--context_pattern {context_pattern} ' \
                  f'--remove_tmp {remove_tmp} ' \
                  f'--tmp_compression {tmp_compression} '
        if max_cov_cutoff is not None:
            command += f'--max_cov_cutoff {max_cov_cutoff}'

        command_dict = {
            'command': command,
            '-pe smp': 1,  # cpu for each command
        }

        command_dict_list.append(command_dict)
    with open(command_file_path, 'w') as f:
        json.dump(command_dict_list, f)

    # Qsub command list
    submitter = Qsubmitter(command_file_path=command_file_path,
                           project_name=project_name,
                           auto_submit=submit)
    return submitter


def _assemble_dataset(cell_meta_df, out_path,
                      region_meta_df, region_name,
                      context_list, cell_region_dir,
                      remove_cell_files=True):
    # qc not happen in this func, assume everything is not missing. only use it inside prepare_dataset()
    h5f = h5py.File(out_path, 'x')  # fail if file exist

    # build region groups, order: region_set, context, cell
    for _region_name, _region_meta in zip(region_name, region_meta_df):
        print(f'Building {_region_name} group...', end=' ')
        region_group = h5f.require_group(_region_name)
        # add region meta
        add_df2group(region_group, 'region_meta', _region_meta)
        # build context group
        for context in context_list:
            print(context, end=' ')
            context_group = region_group.require_group(context)
            # build cell group
            for cell in cell_meta_df.index:
                cell_group = context_group.require_group(cell)
                add_lil_matrix_to_cell_group(cell_group=cell_group,
                                             context=context,
                                             cell_id=cell,
                                             region_name=_region_name,
                                             cell_region_dir=cell_region_dir,
                                             remove_cell_files=remove_cell_files)
        print()
    # add cell meta
    print('Add cell meta into h5mc')
    add_df2group(h5f, 'cell_meta', cell_meta_df)

    h5f.close()
    return


def add_df2group(h5_group, name, df, compression='gzip'):
    # save index into col
    df = df.copy()
    df.reset_index(drop=False, inplace=True)

    # judge dtype
    h5_dtype = []
    for col, col_dt in df.dtypes.iteritems():
        if col_dt == 'O':
            h5_dtype.append((col, h5py.special_dtype(vlen=str)))
        elif col_dt == 'int':
            h5_dtype.append((col, np.int32))
        elif col_dt == 'float':
            h5_dtype.append((col, np.float32))
        else:
            h5_dtype.append((col, col_dt))
    # h5 compound dtype
    h5_dtype = np.dtype(h5_dtype)
    h5_data = [tuple(row) for i, row in df.iterrows()]

    h5_group.require_dataset(name, shape=(df.shape[0],), dtype=h5_dtype, compression=compression)
    h5_group[name][...] = h5_data
    return


def add_lil_matrix_to_cell_group(cell_group, context, cell_id, region_name, cell_region_dir, remove_cell_files):
    file_path = f'{cell_region_dir}/{cell_id}.{region_name}_{context}.count_table.bed.gz'
    df = pd.read_table(file_path, header=None, na_values='.')
    lil_data = lil_matrix(df[[4, 5]].fillna(0).T.values, dtype=uint32)
    mc, cov = lil_data.data
    mc_idx, cov_idx = lil_data.rows
    cell_group.require_dataset('mc', data=[mc, mc_idx], shape=(2, len(mc)), dtype=uint32, compression='gzip')
    cell_group.require_dataset('cov', data=[cov, cov_idx], shape=(2, len(cov)), dtype=uint32, compression='gzip')

    if remove_cell_files:
        run(['rm', '-f', file_path])
    return


def _parse_bed(file_path):
    df = pd.read_table(file_path,
                       header=None,
                       names=['chrom', 'start', 'end', '_id'],
                       index_col='_id')
    return df


def parse_cell_meta_df(cell_meta_path, datasets=None, dataset_col=None, index_col='_id'):
    cell_total_df = pd.read_table(cell_meta_path, header=0, index_col=index_col)
    if dataset_col is not None and datasets is not None:
        if isinstance(datasets, str):
            datasets = [datasets]
        cell_total_df = cell_total_df[cell_total_df[dataset_col].isin(datasets)]
    print('Got %d cells from cell meta table.' % cell_total_df.shape[0])
    return cell_total_df


def prepare_dataset(cell_meta_path, dataset_name, out_dir, genome, cell_id_col='_id',
                    dataset_col=None, allc_path_col='ALLC_path', select_cells=False,
                    region_bed_path=None, region_name=None,
                    context_pattern=None, max_cov_cutoff=None,
                    total_cpu=50, submission_gap=5, qstat_gap=100):
    generation_note = {'assembly_time': cur_time()}

    # prepare refs, change None in to default from config
    genome = genome.upper()
    cell_level_region_set = []
    cell_level_region_name = []
    for _region_set, _region_name in zip(ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP'].split(' '),
                                         ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP_NAME'].split(' ')):
        cell_level_region_set.append(ref_path_config[genome][_region_set])
        cell_level_region_name.append(_region_name)
    if region_bed_path is None:
        region_bed_path = cell_level_region_set
    if region_name is None:
        region_name = cell_level_region_name
    if context_pattern is None:
        context_pattern = ref_path_config['DATA_SELECT']['ALLC_CONTEXT'].split(' ')

    # prepare cell meta df, select dataset or not
    if not select_cells:
        dataset_col = None
        print(f'Using all cells from {cell_meta_path}')
    else:
        if dataset_col is None:
            raise ValueError('dataset_col can not be None, when select_cells is True')
    cell_meta_df = parse_cell_meta_df(cell_meta_path, dataset_name, dataset_col, index_col=cell_id_col)
    # rename required cols for standard: cell id, dataset, allc path
    if dataset_col is None:
        dataset_col = 'dataset'
        cell_meta_df['dataset'] = dataset_name
    cell_meta_df.rename(inplace=True, columns={cell_id_col: '_id',
                                               dataset_col: 'dataset',
                                               allc_path_col: 'ALLC_path'})
    allc_path_col = 'ALLC_path'  # make sure it is the standard name

    # remove cell without allc_path
    cells_no_allc = cell_meta_df[cell_meta_df[allc_path_col].isnull()]
    cell_meta_df = cell_meta_df[~cell_meta_df[allc_path_col].isnull()]
    if cells_no_allc.shape[0] > 0:
        print(f'{cells_no_allc.shape[0]} cell(s) do not have ALLC path, '
              f'removed from cell meta table, see generation note later.')
        generation_note['cell without ALLC'] = cells_no_allc.index.tolist()
        print(f'{cell_meta_df.shape[0]} cells remained.')

    # map to region
    cell_ids = cell_meta_df[allc_path_col].index.tolist()
    allc_files = cell_meta_df[allc_path_col].tolist()
    submitter = batch_map_to_region(cell_ids, allc_files, out_dir, genome, dataset_name,
                                    region_bed_path=region_bed_path, region_name=region_name,
                                    context_pattern=context_pattern, max_cov_cutoff=max_cov_cutoff,
                                    remove_tmp=True, tmp_compression=False,
                                    total_cpu=total_cpu, submission_gap=submission_gap,
                                    qstat_gap=qstat_gap, submit=True)
    # check submitter status before next step
    print(f'{len(submitter.commands)} jobs finished.')

    # make h5 file
    out_path = out_dir + f'/{dataset_name}.h5mc'
    region_meta_df = [_parse_bed(rp) for rp in region_bed_path]
    _assemble_dataset(cell_meta_df=cell_meta_df,
                      out_path=out_path,
                      region_meta_df=region_meta_df,
                      region_name=region_name,
                      context_list=context_pattern,
                      cell_region_dir=out_dir,
                      remove_cell_files=True)

    # write generation note to h5mc file
    # add qsub, ref config and function parameters to h5mc
    with h5py.File(out_path, 'a') as h5f:
        h5f.attrs['qsub_config'] = json.dumps(dict(qsub_config._sections))
        h5f.attrs['reference_config'] = json.dumps(dict(ref_path_config._sections))
        local_var = locals()
        args = inspect.getfullargspec(prepare_dataset).args
        arg_dict = {k: str(v) for k, v in local_var.items() if k in args}
        h5f.attrs['assembly_args'] = json.dumps(arg_dict)
        h5f.attrs['generation_note'] = json.dumps(generation_note)

    return
