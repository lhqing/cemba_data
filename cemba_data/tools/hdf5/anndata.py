"""
BSD 3-Clause License

Copyright (c) 2017 F. Alexander Wolf, P. Angerer, Theis Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import numpy as np
import pandas as pd
from anndata import AnnData
import logging

log = logging.getLogger()


def _get_mean_var(X):
    # - using sklearn.StandardScaler throws an error related to
    #   int to long trafo for very large matrices
    # - using X.multiply is slower
    mean = X.mean(axis=0)
    # scanpy deal with both sparse and full matrix, here only support full
    # if issparse(X):
    #     mean_sq = X.multiply(X).mean(axis=0)
    #     mean = mean.A1
    #     mean_sq = mean_sq.A1
    # else:
    #     mean_sq = np.multiply(X, X).mean(axis=0)
    mean_sq = np.multiply(X, X).mean(axis=0)
    # enforce R convention (unbiased estimator) for variance
    var = (mean_sq - mean**2) * (X.shape[0]/(X.shape[0]-1))
    return mean, var


def highly_variable_feature(
        adata,
        min_disp=None, max_disp=None,
        min_mean=None, max_mean=None,
        n_top_feature=None,
        n_bins=20,
        no_info_cutoff=0.02):
    """
    Adapted from Scanpy, see license above
    """
    log.info('extracting highly variable features')

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`highly_variable_feature` expects an `AnnData` argument, '
            'pass `inplace=False` if you want to return a `np.recarray`.')

    if n_top_feature is not None and not all([min_disp is None, max_disp is None,
                                              min_mean is None, max_mean is None]):
        log.info('If you pass `n_top_feature`, all cutoffs are ignored.')
    if min_disp is None:
        min_disp = 0.5
    if min_mean is None:
        min_mean = 0.
    if max_mean is None:
        max_mean = 10

    # X = np.expm1(adata.X) if flavor == 'seurat' else adata.X
    # not doing any transform for X
    total_X = adata.X
    mean, var = _get_mean_var(total_X)
    # calculate total mean var 1st, filter out those region that have no coverage. (mean -> 1, var -> 0)
    use_var = adata.var.index[(mean - 1) ** 2 + var > no_info_cutoff]
    X = adata[:, use_var].X
    mean, var = _get_mean_var(X)
    # now actually compute the dispersion
    if 0 in mean:
        mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    # raw dispersion is the variance normalized by mean
    dispersion = var / mean
    if 0 in dispersion:
        dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)

    # all of the following quantities are "per-feature" here
    df = pd.DataFrame(index=use_var)
    df['mean'] = mean
    df['dispersion'] = dispersion

    # Seurat's flavor
    df['mean_bin'] = pd.cut(df['mean'], bins=n_bins)
    disp_grouped = df.groupby('mean_bin')['dispersion']
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)
    # retrieve those genes that have nan std, these are the ones where
    # only a single gene fell in the bin and implicitly set them to have
    # a normalized disperion of 1
    one_gene_per_bin = disp_std_bin.isnull()
    gen_indices = np.where(one_gene_per_bin[df['mean_bin']])[0].tolist()
    if len(gen_indices) > 0:
        log.info(
            f'Gene indices {gen_indices} fell into a single bin: their '
            f'normalized dispersion was set to 1.\n'
            f'Decreasing `n_bins` will likely avoid this effect.')
    # Circumvent pandas 0.23 bug. Both sides of the assignment have dtype==float32,
    # but there’s still a dtype error without “.value”.
    disp_std_bin[one_gene_per_bin] = disp_mean_bin[one_gene_per_bin].values
    disp_mean_bin[one_gene_per_bin] = 0
    # actually do the normalization
    df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
                             - disp_mean_bin[df['mean_bin']].values) / disp_std_bin[df['mean_bin']].values
    dispersion_norm = df['dispersion_norm'].values.astype('float32')

    # Select n_top_feature
    if n_top_feature is not None:
        dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
        dispersion_norm[::-1].sort()  # interestingly, np.argpartition is slightly slower
        disp_cut_off = dispersion_norm[n_top_feature - 1]
        gene_subset = np.nan_to_num(df['dispersion_norm'].values) >= disp_cut_off
        log.info(f'the {n_top_feature} top genes correspond to a normalized dispersion cutoff of {disp_cut_off}')
    else:
        max_disp = np.inf if max_disp is None else max_disp
        dispersion_norm[np.isnan(dispersion_norm)] = 0  # similar to Seurat
        gene_subset = np.logical_and.reduce((mean > min_mean, mean < max_mean,
                                             dispersion_norm > min_disp,
                                             dispersion_norm < max_disp))
    df['gene_subset'] = gene_subset
    log.info('    finished')

    # normalize only performed on use_var,
    # add mean and minimum dispersion for those vars not passing no_info_cutoff
    # these region must be filtered out anyway
    other_var = adata[:, ~adata.var.index.isin(use_var)].X
    other_mean, other_var = _get_mean_var(other_var)
    other_df = pd.DataFrame({
        'mean': other_mean,
        'dispersion': df['dispersion'].min(),
        'dispersion_norm': df['dispersion_norm'].min(),
        'gene_subset': False,
    }, index=adata.var.index[~adata.var.index.isin(use_var)])
    df = pd.concat([df, other_df], sort=True).reindex(adata.var.index)
    return df
