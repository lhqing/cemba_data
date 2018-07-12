"""
This file is just for the convenient of myself
"""

from .prepare_dataset import ref_path_config
from .prepare_study import dataset_config
import pandas as pd
import functools


def catch_exception(f):
    @functools.wraps(f)
    def func(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception:
            print('dog is mine, may not work out of lab, plz get your own dog ;)')
    return func


class _Dog:
    def __init__(self):
        self.dataset_config = dataset_config
        self.ref_path_config = ref_path_config
        return

    @catch_exception
    def get_cemba_rs1_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['CEMBA_RS1_METHY']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='region')

    @catch_exception
    def get_human_snmc_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['HUMAN_SNMC']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='dataset')


dog = _Dog()


def _get_cell_metadata_df(cell_meta_path, region_list, dataset_col):
    cell_total_df = pd.read_table(cell_meta_path, header=0, index_col='_id')
    if region_list is not None:
        if isinstance(region_list, str):
            region_list = region_list.split(' ')
        region_list = set([i.upper() for i in region_list])
        cell_total_df = cell_total_df[cell_total_df[dataset_col].isin(region_list)]
    print('Got %d cells from cell meta table.' % cell_total_df.shape[0])
    print('Got %d regions from cell meta table.' % cell_total_df[dataset_col].unique().size)
    return cell_total_df
