import h5py
import numpy as np
import pandas as pd
import json
from ast import literal_eval
from scipy.sparse import lil_matrix, vstack


class Dataset:
    def __init__(self, h5mc_path, mode='r'):
        self.h5f = h5py.File(h5mc_path, mode=mode)
        return

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.h5f.__exit__(exc_type, exc_value, traceback)

    def close(self):
        self.h5f.close()

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

    def iter_cells_data(self, region_name, context, cells, sparse_format='coo'):
        context = context.upper()
        for cell in cells:
            query_path = f'/{region_name}/{context}/{cell}'
            cell_group = self.h5f[query_path]
            cell_mc = _parse_lil_group(cell_group['mc'], sparse_format)
            cell_cov = _parse_lil_group(cell_group['cov'], sparse_format)
            yield cell, cell_mc, cell_cov

    def get_cells_matrix(self, region_name, context, cells=None, sparse_format='csc', to_df=False):
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

