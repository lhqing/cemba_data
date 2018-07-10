import configparser
import os
import datetime
import json
import numpy as np
import pandas as pd
import h5py
import inspect
import argparse
from numpy import uint16, uint32
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
                           auto_submit=submit, total_cpu=total_cpu,
                           submission_gap=submission_gap, qstat_gap=qstat_gap)
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
                _add_lil_matrix_to_cell_group(cell_group=cell_group,
                                              context=context,
                                              cell_id=cell,
                                              region_meta_df=_region_meta,
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


def _add_lil_matrix_to_cell_group(cell_group, context, cell_id, region_meta_df,
                                  region_name, cell_region_dir, remove_cell_files,
                                  data_dtype=uint16, index_dtype=uint32):
    file_path = f'{cell_region_dir}/{cell_id}.{region_name}_{context}.count_table.bed.gz'
    df = pd.read_table(file_path, header=None, na_values='.', index_col='_id',
                       names=['chrom', 'start', 'end', '_id', 'mc', 'cov'])
    # some region maybe dropped during map
    # (e.g. if ALLC don't have chrM, the intersect step will trim chrM regions.)
    # Here make sure every region from region_meta_df is in the index
    df = df.reindex(region_meta_df.index)

    lil_data = lil_matrix(df[['mc', 'cov']].fillna(0).T.values, dtype=uint32)
    mc, cov = lil_data.data
    mc_idx, cov_idx = lil_data.rows

    mc_group = cell_group.require_group('mc')
    mc_group.require_dataset('lil_data', data=mc, shape=(1, len(mc)), dtype=data_dtype, compression='gzip')
    mc_group.require_dataset('lil_index', data=mc_idx, shape=(1, len(mc_idx)), dtype=index_dtype, compression='gzip')
    mc_group.attrs['shape'] = (1, df.shape[0])  # the real shape of sparse matrix

    cov_group = cell_group.require_group('cov')
    cov_group.require_dataset('lil_data', data=cov, shape=(1, len(cov)), dtype=data_dtype, compression='gzip')
    cov_group.require_dataset('lil_index', data=cov_idx, shape=(1, len(cov_idx)), dtype=index_dtype, compression='gzip')
    cov_group.attrs['shape'] = (1, df.shape[0])  # the real shape of sparse matrix
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
                    total_cpu=50, submission_gap=5, qstat_gap=100,
                    remove_cell_files=True, test=False):
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
    cell_meta_df = parse_cell_meta_df(cell_meta_path,
                                      datasets=dataset_name,  # select dataset based on name
                                      dataset_col=dataset_col,
                                      index_col=cell_id_col)
    if test:
        print('Test=True, use the first 10 cells to do test run...')
        if cell_meta_df.shape[0] > 10:
            cell_meta_df = cell_meta_df.iloc[:10, :]
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
                      remove_cell_files=remove_cell_files)

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


def prepare_dataset_register_subparser(subparser):
    parser = subparser.add_parser('prepare-dataset',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Prepare region level methylation dataset for "
                                       "a group of cells defined by a cell metadata table "
                                       "and several reference region sets defined by BED files. \n"
                                       "All data will be saved into a .h5mc file, which is a customized hdf5 file.")
    parser.set_defaults(func=prepare_dataset)

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--cell_meta_path",
        type=str,
        required=True,
        help="Cell metadata table. Format requirement: "
             "1. Must contain a header line. "
             "2. Must have these columns (column name can be specified by corresponding parameters): "
             "   a. cell_id: unique id for each cell; "
             "   b. ALLC_path: ALLC file path for each cell. Cell without ALLC path will be ignored. "
             "3. Optional columns: "
             "   a. dataset: contains the dataset for each cell; "
             "   b. other cell information. "
    )

    parser_req.add_argument(
        "--dataset_name",
        type=str,
        required=True,
        help="Dataset name"
    )

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Output directory"
    )

    parser_req.add_argument(
        "--genome",
        type=str,
        required=True,
        help="Genome used in assembly."
    )

    parser_opt.add_argument(
        "--cell_id_col",
        type=str,
        required=False,
        default='_id',
        help="Name of the cell_id column in cell metadata table."
             " Must specify if different from default."
    )

    parser_opt.add_argument(
        "--dataset_col",
        type=str,
        required=False,
        default='dataset',
        help="Name of the dataset column in cell metadata table. "
             "Must specify if different from default."
    )

    parser_opt.add_argument(
        "--allc_path_col",
        type=str,
        required=False,
        default='ALLC_path',
        help="Name of the ALLC_path column in cell metadata table. "
             "Must specify if different from default."
    )

    parser_opt.add_argument(
        "--select_cells",
        type=str,
        required=False,
        default=False,
        help="If true, will filter cells based on --dataset_name in --dataset_col; \n"
             "If false, will use all cells from the cell metadata table"
    )

    parser_opt.add_argument(
        "--region_bed_path",
        type=str,
        required=False,
        default=None,
        nargs='+',
        help="Space separated region BED file paths for ALLC map-to-region function. "
             "If None, will use default files from config"
    )

    parser_opt.add_argument(
        "--region_name",
        type=str,
        required=False,
        default=None,
        nargs='+',
        help="Space separated region names corresponding to --region_bed_path. "
             "If None, will use default names from config"
    )

    parser_opt.add_argument(
        "--context_pattern",
        type=str,
        required=False,
        default=None,
        nargs='+',
        help="Space separated methylation context pattern, N for ATCG, H for ATC. "
             "If None, will use default context patterns from config"
    )

    parser_opt.add_argument(
        "--max_cov_cutoff",
        type=int,
        required=False,
        default=None,
        help="Maximum cutoff for coverage in each base, "
             "e.g. 2 for single cell data, None for bulk seq."
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        default=None,
        help="Total CPU used in qsub. "
             "If None, will use default from config"
    )

    parser_opt.add_argument(
        "--submission_gap",
        type=int,
        required=False,
        default=None,
        help="Qsub job submission gap. "
             "If None, will use default from config"
    )

    parser_opt.add_argument(
        "--qstat_gap",
        type=int,
        required=False,
        default=None,
        help="Qstat query gap. "
             "If None, will use default from config"
    )

    parser_opt.add_argument(
        "--remove_cell_files",
        type=bool,
        required=False,
        default=True,
        help="Remove cell-region methylation tsv.gz file or not."
    )

    parser_opt.add_argument(
        "--test",
        type=bool,
        required=False,
        default=False,
        help="Test mode, only use the first 10 valid cells to do a test run."
    )

    return
