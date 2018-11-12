"""
dog is public, works everywhere
cat is private, may not work outside my server...

more animals are adding...
"""

from .deprecated_prepare_dataset import ref_path_config
from .deprecated_prepare_study import dataset_config
import pandas as pd
import numpy as np
import collections
import random


def _get_cell_metadata_df(cell_meta_path, region_list, dataset_col, allc_path_col='ALLC_path'):
    cell_total_df = pd.read_table(cell_meta_path, header=0, index_col='_id')
    if region_list is not None:
        if isinstance(region_list, str):
            region_list = region_list.split(' ')
        region_list = set([i.upper() for i in region_list])
        cell_total_df = cell_total_df[cell_total_df[dataset_col].isin(region_list)]
    if cell_total_df[allc_path_col].isnull().sum() != 0:
        print('Cell without ALLC path are removed:',
              cell_total_df[cell_total_df[allc_path_col].isnull()].index.tolist())
        cell_total_df.dropna(subset=[allc_path_col], inplace=True)
    print('Got %d cells from cell meta table.' % cell_total_df.shape[0])
    print('Got %d regions from cell meta table.' % cell_total_df[dataset_col].unique().size)
    return cell_total_df


class _Cat:
    def __init__(self):
        self.dataset_config = dataset_config
        self.ref_path_config = ref_path_config

        self.mouse_gene_gtf = pd.read_table(ref_path_config['MM10']['GENE_FLAT_GTF'],
                                            header=0, index_col='gene_id')
        self.human_gene_gtf = pd.read_table(ref_path_config['HG19']['GENE_FLAT_GTF'],
                                            header=0, index_col='gene_id')

        return

    def get_cemba_rs1_cell_table(self, region_list=None, allc_path_col='ALLC_path'):
        cell_meta_path = self.dataset_config['META_TABLE']['CEMBA_RS1_METHY']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='region',
                                     allc_path_col=allc_path_col)

    def get_human_snmc_cell_table(self, region_list=None, allc_path_col='ALLC_path'):
        cell_meta_path = self.dataset_config['META_TABLE']['HUMAN_SNMC']
        return _get_cell_metadata_df(cell_meta_path, region_list, dataset_col='dataset',
                                     allc_path_col=allc_path_col)

    def get_gene_name(self, gene_id, species=None):
        if len(gene_id) > 6 and gene_id[:7] == 'ENSMUSG':
            return self.mouse_gene_gtf.loc[gene_id, 'gene_name']
        elif len(gene_id) > 4 and gene_id[:4] == 'ENSG':
            return self.human_gene_gtf.loc[gene_id, 'gene_name']
        else:
            # gene name is given, return gene id
            if species is None:
                raise ValueError('Species is None, please specify species when use gene name as input')
            if species.lower() in ['mouse', 'm', 'mm10']:
                sub_gene_df = self.mouse_gene_gtf[self.mouse_gene_gtf['gene_name'] == gene_id]
            elif species.lower() in ['human', 'h', 'hg19', 'hg38']:
                sub_gene_df = self.mouse_gene_gtf[self.human_gene_gtf['gene_name'] == gene_id]
            else:
                raise ValueError('Unknown species', species)
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
