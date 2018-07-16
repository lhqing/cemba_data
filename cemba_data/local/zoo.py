"""
dog is public, works everywhere
cat is private, may not work outside my server...

more animals are adding...
"""

from .prepare_dataset import ref_path_config
from .prepare_study import dataset_config
import pandas as pd
import numpy as np
import collections
import random
import functools


def _catch_exception(f):
    @functools.wraps(f)
    def func(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception:
            print('Cat is mine, may not work out of lab, plz get your own cat ;)')

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

        self.mouse_gene_gtf = pd.read_table(ref_path_config['MM10']['GENE_FLAT_GTF'],
                                            header=0, index_col='gene_id')
        return

    @_catch_exception
    def get_cemba_rs1_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['CEMBA_RS1_METHY']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='region')

    @_catch_exception
    def get_human_snmc_cell_table(self, region_list=None):
        cell_meta_path = self.dataset_config['META_TABLE']['HUMAN_SNMC']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='dataset')

    def get_mouse_gene_name(self, gene_id):
        if len(gene_id) > 6 and gene_id[:7] == 'ENSMUSG':
            return self.mouse_gene_gtf.loc[gene_id, 'gene_name']
        else:
            # gene name is given, return gene id
            sub_gene_df = self.mouse_gene_gtf[self.mouse_gene_gtf['gene_name'] == gene_id]
            if sub_gene_df.shape[0] > 1:
                print(f'Warning, the {gene_id} have {sub_gene_df.shape[0]} matches in gene table.'
                      f' Return the 1st one.')
            return sub_gene_df.index[0]


class _Dog:
    def __init__(self):
        return

    @staticmethod
    def sample_this(obj, n=None, frac=None, axis=None, replace=None, weights=None, seed=1994):
        """
        TODO versatile sample function for any array or matrix like data
        """
        if frac is None and n is None:
            raise ValueError(f"frac and n can't be both None")

        # must be python iterable object
        if not isinstance(obj, collections.Iterable):
            raise TypeError(f"{type(obj)} object is not iterable.")

        if isinstance(obj, (pd.Series, pd.DataFrame)):
            sample = obj.sample(n=n,
                                frac=frac,
                                replace=replace,
                                weights=weights,
                                random_state=seed,
                                axis=axis)
        elif isinstance(obj, np.ndarray):
            if len(obj.shape) == 1:
                # 1D array
                if frac is not None and n is None:
                    n = int(obj.size * frac)
                sample = np.random.choice(obj, size=n, replace=replace, p=weights)
            else:
                # TODO nD array
                print("Out of ability...")
                return
        else:
            if frac is not None and n is None:
                n = int(len(obj) * frac)
            if replace:
                # random choices is with_replacment
                sample = random.choices(list(obj), weights=weights, k=n)
            else:
                sample = random.sample(list(obj), k=n)
        print(f'Random sample size {len(sample)} for object type {type(obj)}, return type {type(sample)}')
        return sample

    @staticmethod
    def chunk_this(obj, bins=None, n=None, axis=None):
        """
        TODO versatile chunk function for any array or matrix like data
        """
        return


dog = _Dog()
cat = _Cat()
