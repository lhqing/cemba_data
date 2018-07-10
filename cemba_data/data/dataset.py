import h5py
import numpy as np
import pandas as pd
import json
from ast import literal_eval
from scipy.sparse import csr_matrix, lil_matrix, vstack


class Dataset:
    def __init__(self, h5mc_path):
        self.h5f_path = h5mc_path
        self.h5f = None
        return

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5f.__exit__(exc_type, exc_value, traceback)

    def close(self):
        self.h5f.close()

    def open(self, mode='r'):
        self.h5f = h5py.File(self.h5f_path, mode=mode)

    def __getitem__(self, item):
        return self.h5f[item]

    @property
    def cell_meta(self):
        meta_df = pd.DataFrame(self.h5f['cell_meta'].value)
        meta_df.set_index('_id', inplace=True)
        return meta_df

    @property
    def assembly_args(self):
        arg_dict = json.loads(self.h5f.attrs['assembly_args'])
        eval_args = {}
        for k, v in arg_dict.items():
            try:
                eval_args[k] = literal_eval(v)
            except (SyntaxError, ValueError):
                eval_args[k] = v
        return eval_args

    def get_region_meta(self, region_name):
        meta_df = pd.DataFrame(self.h5f[region_name]['region_meta'].value)
        meta_df.set_index('_id', inplace=True)
        return meta_df

    # iter_cells_data and get_cells_matrix can be faster...
    def iter_cells_data(self, region_name, context, cells, sparse_format='coo'):
        context = context.upper()
        for cell in cells:
            query_path = f'/{region_name}/{context}/{cell}'
            cell_group = self.h5f[query_path]
            cell_mc = _parse_lil_group(cell_group['mc'], sparse_format)
            cell_cov = _parse_lil_group(cell_group['cov'], sparse_format)
            yield cell, cell_mc, cell_cov

    def get_cells_matrix(self, region_name, context, cells=None, sparse_format='csr', to_df=False):
        region_meta_df = self.get_region_meta(region_name)
        mc_list = []
        cov_list = []
        if cells is None:
            cells = self.cell_meta.index
        for cell, mc, cov in self.iter_cells_data(region_name, context, cells, sparse_format):
            mc_list.append(mc)
            cov_list.append(cov)
        mc_sparse_matrix = vstack(mc_list)  # concatenate by row
        cov_sparse_matrix = vstack(cov_list)
        if to_df:
            print(f'Return pandas DataFrame for {len(cells)} cells')
            mc_df = pd.DataFrame(mc_sparse_matrix.toarray(), columns=region_meta_df.index, index=cells)
            cov_df = pd.DataFrame(cov_sparse_matrix.toarray(), columns=region_meta_df.index, index=cells)
            return mc_df, cov_df
        else:
            column_index = region_meta_df.index
            row_index = cells
            return mc_sparse_matrix, cov_sparse_matrix, column_index, row_index

    # must used func
    def get_mc_rate(self, region_name, context, cov_cutoff, cells=None):
        mc_sparse_matrix, cov_sparse_matrix, column_index, row_index = \
            self.get_cells_matrix(region_name, context, cells)

        mask = cov_sparse_matrix > cov_cutoff
        cov_pass = cov_sparse_matrix.multiply(mask)  # value not pass cutoff will be 0
        mc_pass = mc_sparse_matrix.multiply(mask)  # use same filter as cov
        mc_rate = (mc_pass / cov_pass).astype(np.float32)
        # 0 is different from NA
        mc_rate[mc_rate == 0] = 1e-9  # assign 0 to a small number to be compatible with sparse matrix

        return MethylRate(csr_matrix(mc_rate), column_index, row_index, cov_cutoff)


class MethylRate:
    def __init__(self, mc_rate_csr, col_idx, row_idx, cov_cutoff):
        # input
        self._mc_rate = mc_rate_csr
        self._col_idx = col_idx
        self._row_idx = row_idx
        self._cov_cutoff = cov_cutoff
        # na count
        # turn to csr method
        # self._cell_na_count = np.isnan(mc_rate).sum(axis=1)
        # self._region_na_count = np.isnan(mc_rate).sum(axis=0)
        # mask
        self._cell_mask = None
        self._cell_na_cutoff = None
        self._region_mask = None
        self._region_na_cutoff = None

    def __add__(self, obj):
        pass
        # if not isinstance(obj, MethylRate):
        #    raise TypeError(f'Adding MethylRate objects with {type(obj)} is not allowed.')
        # if self._col_idx != obj._col_idx:
        #    raise ValueError('Adding MethylRate objects with different regions')
        # if self._cov_cutoff != obj._cov_cutoff:
        #    raise ValueError('Adding MethylRate objects with different _cov_cutoff')
        #
        # self._mc_rate = np.concatenate([self._mc_rate, obj._mc_rate], axis=0)
        # self._row_idx = np.concatenate([self._row_idx, obj._row_idx])
        # self._cell_na_count = np.concatenate([self._cell_na_count, obj._cell_na_count])
        # self._region_na_count = self._region_na_count + obj._region_na_count
        ## data changed, must recalculate region mask
        # self._cell_mask = None
        # self._cell_na_cutoff = None
        # self._region_mask = None
        # self._region_na_cutoff = None

    # add a function that allows region concatenate

    @property
    def value(self):
        return self._mc_rate

    def filter_cell(self, cutoff):
        self._cell_na_cutoff = cutoff
        if self._region_mask is None:
            self._cell_mask = _get_na_rate(self._mc_rate, axis=1) < cutoff
        else:
            self._cell_mask = _get_na_rate(self._mc_rate[:, self._region_mask.A1], axis=1) < cutoff

    def filter_region(self, cutoff):
        self._region_na_cutoff = cutoff
        if self._cell_mask is None:
            self._region_mask = _get_na_rate(self._mc_rate, axis=0) < cutoff
        else:
            self._region_mask = _get_na_rate(self._mc_rate[self._cell_mask.A1, :], axis=0) < cutoff


def _get_na_rate(_array, axis):
    pass
    # csr na rate
    #return np.isnan(_array).sum(axis=axis) / _array.shape[axis]


def _parse_lil_group(lil_group, sparse_format):
    data = lil_group['lil_data'].value
    index = lil_group['lil_index'].value
    shape = lil_group.attrs['shape']
    lil = lil_matrix(np.empty(shape=shape), dtype=data.dtype)
    lil.data = data
    lil.rows = index
    if sparse_format == 'coo':
        return lil.tocoo()
    elif sparse_format == 'csr':
        return lil.tocsr()
    elif sparse_format == 'csc':
        return lil.tocsc()
    elif sparse_format == 'lil':
        return lil
    else:
        raise NotImplementedError('sparse_format only support coo, csr, csc, lil.')
