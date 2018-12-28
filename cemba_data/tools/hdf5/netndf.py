import numpy as np
import pandas as pd
import xarray as xr
from anndata import AnnData


def _calculate_posterior_mc_rate(mc_array, cov_array, var_dim='chrom100k',
                                 normalize_per_cell=True, clip_norm_value=10):
    raw_rate = mc_array / cov_array
    cell_rate_mean = raw_rate.mean(dim=var_dim)
    cell_rate_var = raw_rate.var(dim=var_dim)

    # based on beta distribution mean, var
    # a / (a + b) = cell_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = cell_rate_var
    # calculate alpha beta value for each cell
    cell_a = (1 - cell_rate_mean) * (cell_rate_mean ** 2) / cell_rate_var - \
             cell_rate_mean
    cell_b = cell_a * (1 / cell_rate_mean - 1)

    # cell specific posterior rate
    post_rate = (mc_array + cell_a) / (cov_array + cell_a + cell_b)

    if normalize_per_cell:
        # this is normalize by post_rate per cell, just a mean center
        post_rate = post_rate / post_rate.mean(dim='chrom100k')
        if clip_norm_value is not None:
            post_rate.values[np.where(post_rate.values > clip_norm_value)] = clip_norm_value
    return post_rate


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
                    normalize_per_cell=True,
                    clip_norm_value=10,
                    rate_da_suffix='rate'):
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')

        rate = _calculate_posterior_mc_rate(mc_array=da_mc.values,
                                            cov_array=da_cov.values,
                                            normalize_per_cell=normalize_per_cell,
                                            clip_norm_value=clip_norm_value)
        da_rate = xr.DataArray(data=rate,
                               coords=da_mc.coords,
                               dims=da_mc.dims)
        self[da + "_" + rate_da_suffix] = da_rate
        return

    def add_gene_rate(self, dim='gene', da='gene_da',
                      normalize_per_cell=True, clip_norm_value=10,
                      rate_da_suffix='rate'):
        if dim not in self[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        da_mc = self[da].sel(count_type='mc')
        da_cov = self[da].sel(count_type='cov')

        # for gene, we just use normal rate
        rate = da_mc / da_cov

        if normalize_per_cell:
            cell_overall = da_mc.sum(dim='gene') / da_cov.sum(dim='gene')
            rate = rate / cell_overall
            if clip_norm_value is not None:
                rate.values[np.where(rate.values > clip_norm_value)] = clip_norm_value
        self[da + "_" + rate_da_suffix] = rate
        return

    def to_ann(self, da, var_dim, mc_type, obs_dim='cell'):
        index_dict = self[da].indexes
        return AnnData(X=self[da].sel(mc_type=mc_type).values.copy(),
                       obs=pd.DataFrame(index=index_dict[obs_dim]),
                       var=pd.DataFrame(index=index_dict[var_dim]))

    def add_ann_results(self, adata, var_dim, obs_dim='cell'):
        # columns from AnnData.obs and AnnData.var go to da.coords
        # obsm goes to new da with corresponding new dimension
        obs_df = adata.obs
        obs_df.index.name = obs_dim  # make sure coords added with "cell" index
        for col, data in obs_df.iteritems():
            self.coords[col] = data

        var_df = adata.var
        var_df.index.name = var_dim  # make sure coords added with "cell" index
        for col, data in var_df.iteritems():
            self.coords[col] = data

        for obsm_key in adata.obsm_keys():
            coord_name = obsm_key[2:]  # remove 'X_'
            obsm_data = adata.obsm[obsm_key]
            obsm_df = pd.DataFrame(obsm_data,
                                   index=adata.obs_names,
                                   columns=[f'{coord_name}_{i}' for i in range(obsm_data.shape[1])])
            obsm_df.index.name = obs_dim
            obsm_df.columns.name = coord_name
            self[coord_name+'_coord'] = obsm_df

        for varm_key in adata.varm_keys():
            coord_name = varm_key
            varm_data = adata.varm[varm_key]
            varm_df = pd.DataFrame(varm_data,
                                   index=adata.var_names,
                                   columns=[f'{coord_name}_{i}' for i in range(varm_data.shape[1])])
            varm_df.index.name = var_dim
            varm_df.columns.name = coord_name
            self[coord_name + '_coord'] = varm_df
        return

    def add_dataframe_to_coords(self, df, index_dim):
        # add columns to da.coords based on index and index_name
        df.index.name = index_dim
        for col, data in df.iteritems():
            self.coords[col] = data
        return

    def add_dataframe_to_da(self, df, index_dim, col_dim, da_name):
        # add columns to da.coords based on index and index_name
        df.index.name = index_dim
        df.columns.name = col_dim
        self[da_name] = df
        return

    def plotting_df(self):
        return
