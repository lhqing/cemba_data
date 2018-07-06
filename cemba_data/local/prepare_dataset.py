import configparser
import os
import datetime
from .qsub import qsub_config

ref_path_config = configparser.ConfigParser()
ref_path_config.read(os.path.dirname(__file__) + '/config_ref_path.ini')


def cur_time():
    return datetime.datetime.now().strftime("%Y%m%d-%H%M%S")


def batch_map_to_region(data_set_name, cell_ids, allc_files, out_dir, genome,
                        region_bed_path=None, region_name=None, genome_size_path=None,
                        context_pattern=None, max_cov_cutoff=None,
                        remove_tmp=True, tmp_compression=False,
                        total_cpu=50, submission_gap=5,
                        qstat_gap=100):
    # update qsub_config
    qsub_config['RUNNING_DEFAULT']['TOTAL_CPU'] = total_cpu
    qsub_config['RUNNING_DEFAULT']['SUBMISSION_GAP'] = submission_gap
    qsub_config['RUNNING_DEFAULT']['QSTST_GAP'] = qstat_gap
    job_dir = qsub_config['QSUB_DEFAULT']['JOB_DIR']
    project_name = f'map_to_region_{data_set_name}'
    working_job_dir = job_dir + f'/{project_name}'
    try:
        os.mkdir(working_job_dir)
    except OSError:
        pass

    # make region bed info
    genome = genome.upper()
    CELL_LEVEL_REGION_SET = []
    CELL_LEVEL_REGION_NAME = []
    for _region_set, _region_name in zip(ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP'].split(' '),
                                         ref_path_config['DATA_SELECT']['CELL_LEVEL_MAP_NAME'].split(' ')):
        CELL_LEVEL_REGION_SET.append(ref_path_config[genome][_region_set])
        CELL_LEVEL_REGION_NAME.append(ref_path_config[genome][_region_name])

    if region_bed_path is None:
        region_bed_path = CELL_LEVEL_REGION_SET
    if region_name is None:
        region_name = CELL_LEVEL_REGION_NAME
    if genome_size_path is None:
        genome_size_path = ref_path_config[genome]['CHROM_SIZE']
    if context_pattern is None:
        context_pattern =


    # make command json
    command_dict_list = []
    if len(cell_ids) != len(allc_files):
        raise ValueError('Cell id and allc files have different length.')
    for cell, allc in zip(cell_ids, allc_files):
        command_dict = {
            'command': f'cemba_data map-to-region --allc_path {allc} '
                       f'--out_path_prefix {out_dir} '
                       '--region_bed_path {} '
                       '--region_name {} '
                       '--genome_size_path GENOME_SIZE_PATH '
                       '--context_pattern '
                       '--max_cov_cutoff '
                       '--remove_tmp '
                       '--tmp_compression ',
            '-pe smp': 1,  # cpu for each command
        }
    return


def assembl_dataset():
    return


def prepare_dataset():
    # cell df and filter parms (defalut none if df is whole dataset)
    # 1. use batch_map_to_region get all files
    # 2. assemble_dataset to h5 file

    return

