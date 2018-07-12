import h5py
import numpy as np
import pandas as pd
import json
from ast import literal_eval
from scipy.sparse import csr_matrix, lil_matrix, vstack, hstack
from anndata import AnnData


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
        mc_rate[np.isnan(mc_rate)] = 0  # assign 0 to a small number to be compatible with sparse matrix
        return Study(csr_matrix(mc_rate), column_index, row_index, cov_cutoff, context, region_name)


class Study:
    """
    For general analysis, currently I actually use AnnData and scanpy,
    because I don't think its necessary to rebuild the wheels.
    Study can be easily transferred into AnnData by .to_ann()

    The reason I write Study is to host some methods that is specific to methylation data.
    The output file of Study is actually AnnData too,
    which can also be load as a Study using prepare_study.read_from_ann()
    """
    def __init__(self, mc_rate_csr, col_idx, row_idx, cov_cutoff, context, region_name):
        # input
        self._mc_rate = mc_rate_csr
        if isinstance(context, str):
            # only add context suffix when its str (initial study)
            self._col_idx = col_idx + '_' + context
        else:
            self._col_idx = col_idx
        self._row_idx = row_idx

        # because region_append may cause multiple region names, cov cutoffs or contexts
        if isinstance(cov_cutoff, list):
            self._cov_cutoff = cov_cutoff
        else:
            self._cov_cutoff = [cov_cutoff]
        if isinstance(region_name, list):
            self._region_name = region_name
        else:
            self._region_name = [region_name]
        if isinstance(context, list):
            self._context = context
        else:
            self._context = [context]

        # mask
        self._cell_mask = None
        self._cell_na_cutoff = None
        self._region_mask = None
        self._region_na_cutoff = None

    def __add__(self, obj):
        """
        row (cell) concatenate
        :param obj:
        :return:
        """
        if not isinstance(obj, Study):
            raise TypeError(f'Adding MethylRate objects with {type(obj)} is not allowed.')
        if self._col_idx.size != obj._col_idx.size or self._region_name != obj._region_name:
            raise ValueError('Adding MethylRate objects must have same regions')
        if self._cov_cutoff != obj._cov_cutoff:
            raise ValueError('Adding MethylRate objects with different cov_cutoff:',
                             self._cov_cutoff, obj._cov_cutoff)
        if self._context != obj._context:
            raise ValueError('Adding MethylRate objects with different context:',
                             self._context, obj._context)

        new_mc_rate = vstack([self._mc_rate.tocsr(), obj._mc_rate.tocsr()])
        new_row_idx = pd.Series(np.concatenate([self._row_idx, obj._row_idx]))
        # check new row idx, raise if duplicates found
        if new_row_idx.duplicated().sum() != 0:
            raise ValueError('Concatenated cells have duplicate index.')

        return Study(new_mc_rate, self._col_idx, new_row_idx, self._cov_cutoff,
                     self._context, self._region_name)

    def __radd__(self, other):
        return self + other

    # add a function that allows region concatenate
    def region_append(self, obj):
        """
        column (region or context) concatenate
        :param obj: MethylRate object
        :return:
        """
        if not isinstance(obj, Study):
            raise TypeError(f'region_append MethylRate objects with {type(obj)}')
        for x, y in zip(self._row_idx, obj._row_idx):
            if x != y:
                raise ValueError('region_append MethylRate objects must have same cells')
        # region append may allow different cov cutoff
        # if self._cov_cutoff != obj._cov_cutoff:
        #     raise ValueError('region_append MethylRate objects with different cov_cutoff:',
        #                      self._cov_cutoff, obj._cov_cutoff)
        # region append may also allow different context
        # if self._context != obj._context:
        #     raise ValueError('region_append MethylRate objects with different context:',
        #                      self._context, obj._context)

        new_mc_rate = hstack([self._mc_rate.tocsc(), obj._mc_rate.tocsc()])
        new_col_idx = pd.Series(np.concatenate([self._col_idx, obj._col_idx]))

        # check new col idx, raise if duplicates found
        if len(set(new_col_idx)) != len(new_col_idx):
            raise ValueError('Concatenated regions have duplicate names.')

        new_context = self._context + obj._context
        new_region_name = self._region_name + obj._region_name
        new_cov_cutoff = self._cov_cutoff + obj._cov_cutoff
        return Study(new_mc_rate, new_col_idx, self._row_idx, new_cov_cutoff,
                     new_context, new_region_name)

    @property
    def value(self):
        return self._mc_rate

    @property
    def shape(self):
        return self._mc_rate.shape

    def filter_cell(self, max_na_rate):
        self._cell_na_cutoff = max_na_rate
        if self._region_mask is None:
            self._cell_mask = _get_sparse_na_mask(self._mc_rate, axis=1, cutoff=max_na_rate)
        else:
            self._cell_mask = _get_sparse_na_mask(self._mc_rate[:, np.ravel(self._region_mask)],
                                                  axis=1, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._cell_mask) / len(self._cell_mask) * 100)
        print(f'{sum(self._cell_mask)} cells ({na_rate}%) remained after this cell filter')

    def filter_region(self, max_na_rate):
        self._region_na_cutoff = max_na_rate
        if self._cell_mask is None:
            self._region_mask = _get_sparse_na_mask(self._mc_rate, axis=0, cutoff=max_na_rate)
        else:
            self._region_mask = _get_sparse_na_mask(self._mc_rate[np.ravel(self._cell_mask), :],
                                                    axis=0, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._region_mask) / len(self._region_mask) * 100)
        print(f'{sum(self._region_mask)} regions ({na_rate}%) remained after this region filter')

    def reset_filter(self):
        self._region_mask = None
        self._region_na_cutoff = None
        self._cell_mask = None
        self._cell_na_cutoff = None

    def to_ann(self):
        rows = {'row_names': self._row_idx}
        cols = {'col_names': self._col_idx}
        uns = {}
        if self._region_na_cutoff is not None:
            uns['region_na_cutoff'] = self._region_na_cutoff
            cols['region_mask'] = self._region_mask
        if self._cell_na_cutoff is not None:
            uns['cell_na_cutoff'] = self._cell_na_cutoff
            rows['cell_mask'] = self._cell_mask
        return AnnData(self._mc_rate, rows, cols, uns=uns)

    def save(self):
        return


class StudyExp:
    """
    For general analysis, currently I actually use AnnData and scanpy,
    because I don't think its necessary to rebuild the wheels.
    Study can be easily transferred into AnnData by .to_ann()

    The reason I write Study is to host some methods that are specific to methylation data.
    The output file of Study is actually AnnData too,
    which can also be load as a Study using prepare_study.read_from_ann()
    """
    def __init__(self, mc_rate_csr, row_dict, col_dict, uns_dict):
        # input
        self._mc_rate = mc_rate_csr
        if isinstance(context, str):
            # only add context suffix when its str (initial study)
            self._col_idx = col_idx + '_' + context
        else:
            self._col_idx = col_idx
        self._row_idx = row_idx

        # because region_append may cause multiple region names, cov cutoffs or contexts
        if isinstance(cov_cutoff, list):
            self._cov_cutoff = cov_cutoff
        else:
            self._cov_cutoff = [cov_cutoff]
        if isinstance(region_name, list):
            self._region_name = region_name
        else:
            self._region_name = [region_name]
        if isinstance(context, list):
            self._context = context
        else:
            self._context = [context]

        # mask
        self._cell_mask = None
        self._cell_na_cutoff = None
        self._region_mask = None
        self._region_na_cutoff = None

    def __add__(self, obj):
        """
        row (cell) concatenate
        :param obj:
        :return:
        """
        if not isinstance(obj, Study):
            raise TypeError(f'Adding MethylRate objects with {type(obj)} is not allowed.')
        if self._col_idx.size != obj._col_idx.size or self._region_name != obj._region_name:
            raise ValueError('Adding MethylRate objects must have same regions')
        if self._cov_cutoff != obj._cov_cutoff:
            raise ValueError('Adding MethylRate objects with different cov_cutoff:',
                             self._cov_cutoff, obj._cov_cutoff)
        if self._context != obj._context:
            raise ValueError('Adding MethylRate objects with different context:',
                             self._context, obj._context)

        new_mc_rate = vstack([self._mc_rate.tocsr(), obj._mc_rate.tocsr()])
        new_row_idx = pd.Series(np.concatenate([self._row_idx, obj._row_idx]))
        # check new row idx, raise if duplicates found
        if new_row_idx.duplicated().sum() != 0:
            raise ValueError('Concatenated cells have duplicate index.')

        return Study(new_mc_rate, self._col_idx, new_row_idx, self._cov_cutoff,
                     self._context, self._region_name)

    def __radd__(self, other):
        return self + other

    # add a function that allows region concatenate
    def region_append(self, obj):
        """
        column (region or context) concatenate
        :param obj: MethylRate object
        :return:
        """
        if not isinstance(obj, Study):
            raise TypeError(f'region_append MethylRate objects with {type(obj)}')
        for x, y in zip(self._row_idx, obj._row_idx):
            if x != y:
                raise ValueError('region_append MethylRate objects must have same cells')
        # region append may allow different cov cutoff
        # if self._cov_cutoff != obj._cov_cutoff:
        #     raise ValueError('region_append MethylRate objects with different cov_cutoff:',
        #                      self._cov_cutoff, obj._cov_cutoff)
        # region append may also allow different context
        # if self._context != obj._context:
        #     raise ValueError('region_append MethylRate objects with different context:',
        #                      self._context, obj._context)

        new_mc_rate = hstack([self._mc_rate.tocsc(), obj._mc_rate.tocsc()])
        new_col_idx = pd.Series(np.concatenate([self._col_idx, obj._col_idx]))

        # check new col idx, raise if duplicates found
        if len(set(new_col_idx)) != len(new_col_idx):
            raise ValueError('Concatenated regions have duplicate names.')

        new_context = self._context + obj._context
        new_region_name = self._region_name + obj._region_name
        new_cov_cutoff = self._cov_cutoff + obj._cov_cutoff
        return Study(new_mc_rate, new_col_idx, self._row_idx, new_cov_cutoff,
                     new_context, new_region_name)

    @property
    def value(self):
        return self._mc_rate

    @property
    def shape(self):
        return self._mc_rate.shape

    def filter_cell(self, max_na_rate):
        self._cell_na_cutoff = max_na_rate
        if self._region_mask is None:
            self._cell_mask = _get_sparse_na_mask(self._mc_rate, axis=1, cutoff=max_na_rate)
        else:
            self._cell_mask = _get_sparse_na_mask(self._mc_rate[:, np.ravel(self._region_mask)],
                                                  axis=1, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._cell_mask) / len(self._cell_mask) * 100)
        print(f'{sum(self._cell_mask)} cells ({na_rate}%) remained after this cell filter')

    def filter_region(self, max_na_rate):
        self._region_na_cutoff = max_na_rate
        if self._cell_mask is None:
            self._region_mask = _get_sparse_na_mask(self._mc_rate, axis=0, cutoff=max_na_rate)
        else:
            self._region_mask = _get_sparse_na_mask(self._mc_rate[np.ravel(self._cell_mask), :],
                                                    axis=0, cutoff=max_na_rate)
        na_rate = '%.1f' % (sum(self._region_mask) / len(self._region_mask) * 100)
        print(f'{sum(self._region_mask)} regions ({na_rate}%) remained after this region filter')

    def reset_filter(self):
        self._region_mask = None
        self._region_na_cutoff = None
        self._cell_mask = None
        self._cell_na_cutoff = None

    def to_ann(self):
        rows = {'row_names': self._row_idx}
        cols = {'col_names': self._col_idx}
        uns = {}
        if self._region_na_cutoff is not None:
            uns['region_na_cutoff'] = self._region_na_cutoff
            cols['region_mask'] = self._region_mask
        if self._cell_na_cutoff is not None:
            uns['cell_na_cutoff'] = self._cell_na_cutoff
            rows['cell_mask'] = self._cell_mask
        return AnnData(self._mc_rate, rows, cols, uns=uns)

    def save(self):
        return


def _get_sparse_na_mask(sparse_arr, axis, cutoff):
    na_mask = np.isnan(sparse_arr.toarray()).sum(axis=axis) / sparse_arr.shape[axis] < cutoff
    return na_mask


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

