import numpy as np


class MCDS:
    def __init__(self, ds):
        self.dataset = ds
        return

    def filter_cell_cov(self, dim, da, mc_type,
                        min_cov=10000, max_cov=None, inplace=False):
        if dim not in self.dataset[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        cell_sum = self \
            .dataset \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .sum(dim)[da]
        cell_max = cell_sum < max_cov
        cell_min = cell_sum > min_cov
        cell_mask = np.all(np.vstack([cell_max.values,
                                      cell_min.values]),
                           axis=0)
        select_index = self.dataset.get_index('cell')[cell_mask]
        filtered_ds = self.dataset.loc[dict(cell=select_index)]
        if inplace:
            self.dataset = filtered_ds
            return MCDS(filtered_ds)

    def filter_region_cov(self, dim, da, mc_type,
                          min_cov=10000, max_cov=None, inplace=False):
        if dim not in self.dataset[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        region_sum = self \
            .dataset \
            .sel(count_type='cov', mc_type=mc_type) \
            .squeeze() \
            .sum('cell')[da]
        region_max = region_sum < max_cov
        region_min = region_sum > min_cov
        region_mask = np.all(np.vstack([region_max.values,
                                        region_min.values]),
                             axis=0)
        select_index = self.dataset.get_index(dim)[region_mask]
        filtered_ds = self.dataset.loc[dict(cell=select_index)]
        if inplace:
            self.dataset = filtered_ds
            return MCDS(filtered_ds)

    def get_mc_rate(self, dim, da,
                    normalize_cell_overall=True,
                    save_rate_da=True,
                    rate_da_suffix='rate'):
        if dim not in self.dataset[da].dims:
            raise ValueError(f'{dim} is not a dimension of {da}')
        da_mc = self.dataset[da].sel(count_type='mc')
        da_cov = self.dataset[da].sel(count_type='cov')
        da_mc_rate = (da_mc / da_cov)
        if normalize_cell_overall:
            cell_overall_rate = da_mc.sum(dim) / da_cov.sum(dim)
            da_mc_rate /= cell_overall_rate
        if save_rate_da:
            self.dataset[da+"_"+rate_da_suffix] = da_mc_rate
        else:
            return da_mc_rate

    def to_ann(self):
        return

    def gather_results_from_ann(self):
        return

    def plotting_df(self):
        return
