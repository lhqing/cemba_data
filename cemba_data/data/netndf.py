import numpy as np
import pandas as pd
import xarray as xr
from anndata import AnnData


class MCDS(xr.Dataset):
    def __init__(self, dataset):
        super().__init__(data_vars=dataset.data_vars, coords=dataset.coords,
                         attrs=dataset.attrs, compat='broadcast_equals')
        return

    def filter_cell_cov(self, dim, da, mc_type,
                        min_cov=10000, max_cov=None):
        """
        Filter cell by total cov for certain mc_type along certain dimension in certain dataarray.

        Parameters
        ----------
        dim
            region dimension to sum
        da
            dataarray to do calculation
        mc_type
            mc_type to sum
        min_cov
            minimum cov sum, suggest to plot distribution first.
        max_cov
            maximum cov sum, suggest ot plot distribution first.

        Returns
        -------

        """
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        cell_sum = self[da] \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .sum(dim)
        if max_cov is None:
            max_cov = np.inf
        cell_max = cell_sum < max_cov
        cell_min = cell_sum > min_cov
        cell_mask = np.all(np.vstack([cell_max.values,
                                      cell_min.values]),
                           axis=0)
        select_index = self.get_index('cell')[cell_mask]
        return self.loc[dict(cell=select_index)]

    def filter_region_cov(self, dim, da, mc_type,
                          min_cov=None, max_cov=None):
        """
        Filter cell by total cov for certain mc_type along certain dimension in certain dataarray.

        Parameters
        ----------
        dim
            region dimension to filter,
            Note when filtering region, sum is always performed on cell.
        da
            dataarray to do calculation
        mc_type
            mc_type to sum
        min_cov
            minimum cov sum, suggest to plot distribution first.
        max_cov
            maximum cov sum, suggest ot plot distribution first.

        Returns
        -------

        """
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        region_sum = self[da] \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .sum('cell')
        if max_cov is None:
            max_cov = np.inf
        region_max = region_sum < max_cov
        region_min = region_sum > min_cov
        region_mask = np.all(np.vstack([region_max.values,
                                        region_min.values]),
                             axis=0)
        select_index = self.get_index(dim)[region_mask]
        filtered_ds = self.loc[{dim: select_index}]
        return MCDS(filtered_ds)

    def add_mc_rate(self, dim, da,
                    normalize_cell_overall=True,
                    cell_overall_rate=None,
                    rate_da_suffix='rate'):
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')
        da_mc_rate = (da_mc / da_cov)
        if normalize_cell_overall:
            if cell_overall_rate is None:
                cell_overall_rate = da_mc.sum(dim) / da_cov.sum(dim)
            da_mc_rate /= cell_overall_rate
        self[da + "_" + rate_da_suffix] = da_mc_rate

    def to_ann(self, da, var_dim, mc_type, obs_dim='cell'):
        index_dict = self[da].indexes
        return AnnData(X=self[da].sel(mc_type=mc_type).values.copy(),
                       obs=pd.DataFrame(index=index_dict[obs_dim]),
                       var=pd.DataFrame(index=index_dict[var_dim]))

    def from_ann(self):
        return

    def plotting_df(self):
        return
