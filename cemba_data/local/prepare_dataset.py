import configparser
import os
import datetime
import json
from .qsub import Qsubmitter, qsub_config

ref_path_config = configparser.ConfigParser()
ref_path_config.read(os.path.dirname(__file__) + '/config_ref_path.ini')


def cur_time():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def batch_map_to_region(cell_ids, allc_files, out_dir, genome, data_set_name=None,
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
    if data_set_name is None:  # if dataset name not provided, use time
        data_set_name = cur_time()
    project_name = f'map_to_region_{data_set_name}'
    working_job_dir = job_dir + f'/{project_name}'
    command_file_path = working_job_dir + '/command_list.json'
    try:
        os.mkdir(working_job_dir)
    except OSError:
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
    Qsubmitter(command_file_path=command_file_path,
               project_name=project_name,
               auto_submit=submit)
    return Qsubmitter


def assemble_dataset(cell_meta, region_meta, region_names, cell_region_dir, remove_cell_files=True):
    return


def prepare_dataset():
    # cell df and filter parms (defalut none if df is whole dataset)
    # 1. use batch_map_to_region get all files
    # 2. assemble_dataset to h5 file
    # batch_map_to_region()
    # assemble_dataset()
    return
