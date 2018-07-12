import configparser
import os
import json
import pandas as pd
from ..data.hdf5 import Dataset


dataset_config = configparser.ConfigParser()
dataset_config.read(os.path.dirname(__file__) + '/config_dataset_path.ini')
CURRENT_PROJECT = list(dataset_config.keys())
CURRENT_PROJECT.remove('DEFAULT')


def print_config():
    print(json.dumps(dataset_config._sections, indent=4))
    return


def prepare_study(project_name, study_name, cell_list, region,
                  region_context, coverage_cutoff, out_dir,
                  cell_id_dict=None, dataset_path_dict=None):

    # parse cell_list, select dataset
    if cell_id_dict is None:
        # cell_id form: {dataset_name}_certain-id
        # cell_id_dict = {dataset: [cell_id]}
        dataset_series = pd.Series({cell: cell.split('_')[0] for cell in cell_list})
        cell_id_dict = {i: v.index.tolist() for i, v in dataset_series.groupby(dataset_series)}
    else:
        print('cell_id_dict is passed, will ignore cell_list.')

    if dataset_path_dict is not None:
        for k in cell_id_dict.keys():
            if k not in dataset_path_dict:
                if project_name in dataset_config:
                    try:
                        dataset_path_dict[k] = dataset_config[project_name][k]
                    except KeyError:
                        raise KeyError(f'Dataset {k} not found in neither dataset_path_dict nor dataset_config')
                else:
                    raise KeyError(f'Dataset {k} not found in neither dataset_path_dict nor dataset_config')
    else:
        if project_name not in dataset_config:
            raise KeyError(f'Project name {project_name} not in dataset_config. '
                           f'Available projects are: {CURRENT_PROJECT}.')
        try:
            dataset_path_dict = {dataset_config[project_name][k] for k in cell_id_dict.keys()}
        except KeyError:
            raise KeyError(f'Dataset {k} not found in dataset_config')


    return


def _get_study_from_one_dataset(dataset_path, cell_id_list, region, region_context, coverage_cutoff):
    ds = Dataset(dataset_path, 'r')
    ds.get_cells_matrix()
    yield