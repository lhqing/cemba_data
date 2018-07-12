"""
dog is public, works everywhere
cat is private, may not work outside my server...

more animals are adding...
"""

from .prepare_dataset import ref_path_config
from .prepare_study import dataset_config
from multiprocessing import Pool
import pandas as pd
import functools


def _catch_exception(f):
    @functools.wraps(f)
    def func(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception:
            print('dog is mine, may not work out of lab, plz get your own dog ;)')
    return func


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


class _Cat:
    def __init__(self):
        self.dataset_config = dataset_config
        self.ref_path_config = ref_path_config
        return

    @_catch_exception
    def get_cemba_rs1_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['CEMBA_RS1_METHY']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='region')

    @_catch_exception
    def get_human_snmc_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['HUMAN_SNMC']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='dataset')


class _Dog:
    def __init__(self):
        return

    def sample_this(self, obj, rate=None, n=None, axis=None, seed=1994):
        """
        TODO versatile sample function for any array or matrix like data
        """
        return

    def chunk_this(self, obj, bins=None, n=None, axis=None):
        """
        TODO versatile chunk function for any array or matrix like data
        """
        return

    def parallel_this(self, f, cpu):
        @functools.wraps(f)
        def func(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except Exception:
                print('Dog can not parallel this... Run again see what is wrong...')

        return func


dog = _Dog()
cat = _Cat()
