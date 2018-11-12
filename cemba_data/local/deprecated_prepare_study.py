import configparser
import os
import json
import pandas as pd
from functools import reduce
from operator import add
from ..data.hdf5 import Dataset


dataset_config = configparser.ConfigParser()
dataset_config.read(os.path.dirname(__file__) + '/deprecated_config_dataset_path.ini')
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
            dataset_path_dict = {k: dataset_config[project_name][k] for k in cell_id_dict.keys()}
        except KeyError:
            print(project_name)
            print(cell_id_dict.keys())
            raise KeyError(f'Dataset not found in dataset_config')

    if isinstance(region, str):
        region = [region]
    if isinstance(region_context, str):
        region_context = [region_context]
    if isinstance(coverage_cutoff, int):
        coverage_cutoff = [coverage_cutoff]

    if len(region_context) != len(region) or len(region_context) != len(coverage_cutoff):
        raise ValueError('Region context and region and coverage cutoff must have same length.')

    study = None
    for _study in _get_study_from_datasets_dif_col(dataset_path_dict,
                                                   cell_id_dict,
                                                   region,
                                                   region_context,
                                                   coverage_cutoff):
        if study is None:
            study = _study
        else:
            study = study.region_append(_study)

    # TODO: direct save to the out path
    return study


def read_from_ann(ann_path):
    # TODO: this func should only take ann_path as input, try to get all other info from the hdf5
    return


def prepare_study_register_subparser():
    # TODO: parser for prepare study
    return


def _get_study_from_datasets_same_col(dataset_path_dict, cell_id_dict, region, region_context, coverage_cutoff):
    """
    generator of study with same columns, used for add
    """
    for dataset, cell_id_list in cell_id_dict.items():
        ds_path = dataset_path_dict[dataset]
        ds = Dataset(ds_path, 'r')
        study = ds.get_mc_rate(region_name=region,
                               context=region_context,
                               cov_cutoff=coverage_cutoff,
                               cells=cell_id_list)
        yield study


def _get_study_from_datasets_dif_col(dataset_path_dict, cell_id_dict, region, region_context, coverage_cutoff):
    """
    generator of study with different columns, used for region_append
    """
    for _region, _region_context, _coverage_cutoff in zip(region, region_context, coverage_cutoff):
        yield reduce(add, [study for study in _get_study_from_datasets_same_col(
            dataset_path_dict, cell_id_dict, _region, _region_context, _coverage_cutoff
        )])
