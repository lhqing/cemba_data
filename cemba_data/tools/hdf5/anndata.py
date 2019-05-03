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
from ..utilities import get_mean_var
from scipy import stats
from statsmodels.stats.multitest import multipletests
import logging
from math import sqrt, floor
from scipy.sparse import issparse

log = logging.getLogger()


# TODO make a MCAD class that overwrite AnnData function and use OOP in clustering
# put recipe in this class
class MCAD(AnnData):
    def __init__(self, adata):
        super().__init__(X=adata.X,
                         obs=adata.obs, var=adata.var, uns=adata.uns,
                         obsm=adata.obsm, varm=adata.varm,
                         layers=adata.layers, raw=adata.raw,
                         dtype=adata.dtype, shape=adata.shape,
                         filename=adata.filename, filemode=adata.filemode)
        return


def highly_variable_methylation_feature(
        adata, feature_mean_cov,
        min_disp=0.5, max_disp=None,
        min_mean=0, max_mean=5,
        n_top_feature=None, bin_min_features=5,
        mean_binsize=0.05, cov_binsize=100):
    """
    Adapted from Scanpy, see license above
    The main difference is that, this function normalize dispersion based on both mean and cov bins.
    """
    # RNA is only scaled by mean, but methylation is scaled by both mean and cov
    log.info('extracting highly variable features')

    if not isinstance(adata, AnnData):
        raise ValueError(
            '`highly_variable_methylation_feature` expects an `AnnData` argument')

    if n_top_feature is not None:
        log.info('If you pass `n_top_feature`, all cutoffs are ignored.')

    # warning for extremely low cov
    low_cov_portion = (feature_mean_cov < 30).sum() / feature_mean_cov.size
    if low_cov_portion > 0.2:
        log.warning(f'{int(low_cov_portion * 100)}% feature with < 10 mean cov, '
                    f'consider filter by cov before find highly variable feature. '
                    f'Otherwise some low coverage feature may be elevated after normalization.')

    X = adata.X
    cov = feature_mean_cov

    mean, var = get_mean_var(X)
    # now actually compute the dispersion
    if 0 in mean:
        mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    # raw dispersion is the variance normalized by mean
    dispersion = var / mean
    if 0 in dispersion:
        dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)

    # all of the following quantities are "per-feature" here
    df = pd.DataFrame(index=adata.var_names)
    df['mean'] = mean
    df['dispersion'] = dispersion
    df['cov'] = cov

    # instead of n_bins, use bin_size, because cov and mc are in different scale
    df['mean_bin'] = (df['mean'] / mean_binsize).astype(int)
    df['cov_bin'] = (df['cov'] / cov_binsize).astype(int)

    # save bin_count df, gather bins with more than bin_min_features features
    bin_count = df.groupby(['mean_bin', 'cov_bin']) \
        .apply(lambda i: i.shape[0]) \
        .reset_index() \
        .sort_values(0, ascending=False)
    bin_count.head()
    bin_more_than = bin_count[bin_count[0] > bin_min_features]
    if bin_more_than.shape[0] == 0:
        raise ValueError(f'No bin have more than {bin_min_features} features, uss larger bin size.')

    # for those bin have too less features, merge them with closest bin in manhattan distance
    # this usually don't cause much difference (a few hundred features), but the scatter plot will look more nature
    index_map = {}
    for _, (mean_id, cov_id, count) in bin_count.iterrows():
        if count > 1:
            index_map[(mean_id, cov_id)] = (mean_id, cov_id)
        manhattan_dist = (bin_more_than['mean_bin'] - mean_id).abs() + (bin_more_than['cov_bin'] - cov_id).abs()
        closest_more_than = manhattan_dist.sort_values().index[0]
        closest = bin_more_than.loc[closest_more_than]
        index_map[(mean_id, cov_id)] = tuple(closest.tolist()[:2])
    # apply index_map to original df
    raw_bin = df[['mean_bin', 'cov_bin']].set_index(['mean_bin', 'cov_bin'])
    raw_bin['use_mean'] = pd.Series(index_map).apply(lambda i: i[0])
    raw_bin['use_cov'] = pd.Series(index_map).apply(lambda i: i[1])
    df['mean_bin'] = raw_bin['use_mean'].values
    df['cov_bin'] = raw_bin['use_cov'].values

    # calculate bin mean and std, now disp_std_bin shouldn't have NAs
    disp_grouped = df.groupby(['mean_bin', 'cov_bin'])['dispersion']
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)

    # actually do the normalization
    _mean_norm = disp_mean_bin.loc[list(zip(df['mean_bin'], df['cov_bin']))]
    _std_norm = disp_std_bin.loc[list(zip(df['mean_bin'], df['cov_bin']))]
    df['dispersion_norm'] = (df['dispersion'].values  # use values here as index differs
                             - _mean_norm.values) / _std_norm.values
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
    return df


def _select_groups(adata, groups_order_subset='all', key='groups'):
    """Get subset of groups in adata.obs[key].
    """
    groups_order = adata.obs[key].cat.categories
    if key + '_masks' in adata.uns:
        groups_masks = adata.uns[key + '_masks']
    else:
        groups_masks = np.zeros((len(adata.obs[key].cat.categories),
                                 adata.obs[key].values.size), dtype=bool)
        for iname, name in enumerate(adata.obs[key].cat.categories):
            # if the name is not found, fallback to index retrieval
            if adata.obs[key].cat.categories[iname] in adata.obs[key].values:
                mask = adata.obs[key].cat.categories[iname] == adata.obs[key].values
            else:
                mask = str(iname) == adata.obs[key].values
            groups_masks[iname] = mask
    if groups_order_subset != 'all':
        groups_ids = []
        for name in groups_order_subset:
            groups_ids.append(
                np.where(adata.obs[key].cat.categories.values == name)[0][0])
        if len(groups_ids) == 0:
            # fallback to index retrieval
            groups_ids = np.where(
                np.in1d(np.arange(len(adata.obs[key].cat.categories)).astype(str),
                        np.array(groups_order_subset)))[0]
        if len(groups_ids) == 0:
            log.info(np.array(groups_order_subset),
                     'invalid! specify valid groups_order (or indices) one of',
                     adata.obs[key].cat.categories)
            from sys import exit
            exit(0)
        groups_masks = groups_masks[groups_ids]
        groups_order_subset = adata.obs[key].cat.categories[groups_ids].values
    else:
        groups_order_subset = groups_order.values
    return groups_order_subset, groups_masks


def rank_features_groups(
        adata,
        groupby,
        use_raw=True,
        groups='all',
        reference='rest',
        n_genes=100,
        rankby_abs=False,
        copy=False,
        method='t-test_overestim_var',
        corr_method='benjamini-hochberg'):
    """Rank genes for characterizing groups.
    Adapted from scanpy, the mainly modification is that, in scRNA, we try to find HEG,
    but in methylation, we try to find hypo-methylation
    """
    # https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html (marker finding part)

    log.info('ranking genes')
    avail_methods = {'t-test', 't-test_overestim_var', 'wilcoxon'}
    if method not in avail_methods:
        raise ValueError('Method must be one of {}.'.format(avail_methods))

    avail_corr = {'benjamini-hochberg', 'bonferroni'}
    if corr_method not in avail_corr:
        raise ValueError('Correction method must be one of {}.'.format(avail_corr))

    adata = adata.copy() if copy else adata
    # sanitize anndata
    # adata._sanitize()

    # for clarity, rename variable
    groups_order = groups if isinstance(groups, str) else list(groups)
    if isinstance(groups_order, list) and isinstance(groups_order[0], int):
        groups_order = [str(n) for n in groups_order]
    if reference != 'rest' and reference not in set(groups_order):
        groups_order += [reference]
    if (reference != 'rest'
            and reference not in set(adata.obs[groupby].cat.categories)):
        raise ValueError('reference = {} needs to be one of groupby = {}.'
                         .format(reference,
                                 adata.obs[groupby].cat.categories.tolist()))

    groups_order, groups_masks = _select_groups(
        adata, groups_order, groupby)

    key_added = 'rank_genes_groups'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = {
        'groupby': groupby,
        'reference': reference,
        'method': method,
        'use_raw': use_raw,
        'corr_method': corr_method,
    }

    # adata_comp mocks an AnnData object if use_raw is True
    # otherwise it's just the AnnData object
    adata_comp = adata
    if adata.raw is not None and use_raw:
        adata_comp = adata.raw
    X = adata_comp.X

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    if n_genes_user > X.shape[1]:
        n_genes_user = X.shape[1]
    # in the following, n_genes is simply another name for the total number of genes
    n_genes = X.shape[1]

    n_groups = groups_masks.shape[0]
    ns = np.zeros(n_groups, dtype=int)
    for imask, mask in enumerate(groups_masks):
        ns[imask] = np.where(mask)[0].size
    if reference != 'rest':
        ireference = np.where(groups_order == reference)[0][0]
    else:
        ireference = None
    reference_indices = np.arange(adata_comp.n_vars, dtype=int)

    rankings_gene_scores = []
    rankings_gene_names = []
    rankings_gene_logfoldchanges = []
    rankings_gene_pvals = []
    rankings_gene_pvals_adj = []

    if method in {'t-test', 't-test_overestim_var'}:
        # loop over all masks and compute means, variances and sample numbers
        means = np.zeros((n_groups, n_genes))
        vars = np.zeros((n_groups, n_genes))
        for imask, mask in enumerate(groups_masks):
            means[imask], vars[imask] = get_mean_var(X[mask])
        # test each either against the union of all other groups or against a
        # specific group
        for igroup in range(n_groups):
            if reference == 'rest':
                mask_rest = ~groups_masks[igroup]
            else:
                if igroup == ireference:
                    continue
                else:
                    mask_rest = groups_masks[ireference]
            mean_rest, var_rest = get_mean_var(X[mask_rest])
            ns_group = ns[igroup]  # number of observations in group

            # this is the only place t-test different from t-test_overestim_var
            if method == 't-test':
                ns_rest = np.where(mask_rest)[0].size
            elif method == 't-test_overestim_var':
                ns_rest = ns[igroup]  # hack for overestimating the variance for small groups
            else:
                raise ValueError('Method does not exist.')

            # if ns_rest became smaller (t-test_overestim_var),
            # the denominator is larger, t is smaller
            denominator = np.sqrt(vars[igroup] / ns_group + var_rest / ns_rest)
            denominator[np.flatnonzero(denominator == 0)] = np.nan  # set denominator == 0 to nan

            # for methylation, we reverse the score, this doesn't impact p-value, but impact score sorting.
            # because we are looking for hypo-methylation
            scores = (mean_rest - means[igroup]) / denominator
            # for RNA, the original scanpy score is here:
            # scores = (means[igroup] - mean_rest) / denominator #Welch t-test

            mean_rest[mean_rest == 0] = 1e-9  # set 0s to small value
            foldchanges = (means[igroup] + 1e-9) / mean_rest
            scores[np.isnan(scores)] = 0

            # Get p-values
            # https://en.wikipedia.org/wiki/Welch%27s_t-test
            denominator_dof = (np.square(vars[igroup]) / (np.square(ns_group) * (ns_group - 1))) + (
                (np.square(var_rest) / (np.square(ns_rest) * (ns_rest - 1))))
            denominator_dof[np.flatnonzero(denominator_dof == 0)] = np.nan
            dof = np.square(
                vars[igroup] / ns_group + var_rest / ns_rest) / denominator_dof  # dof calculation for Welch t-test
            dof[np.isnan(dof)] = 0
            pvals = stats.t.sf(abs(scores), dof) * 2  # *2 because of two-tailed t-test

            if corr_method == 'benjamini-hochberg':
                pvals[np.isnan(pvals)] = 1  # set Nan values to 1 to properly convert using Benhjamini Hochberg
                _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
            elif corr_method == 'bonferroni':
                pvals_adj = pvals * n_genes
            else:
                raise ValueError('Method does not exist.')

            # this rank rule is same as scanpy, but I modified the score calculation above
            scores_sort = np.abs(scores) if rankby_abs else scores
            partition = np.argpartition(scores_sort, -n_genes_user)[-n_genes_user:]
            partial_indices = np.argsort(scores_sort[partition])[::-1]
            global_indices = reference_indices[partition][partial_indices]
            rankings_gene_scores.append(scores[global_indices])
            rankings_gene_logfoldchanges.append(np.log2(np.abs(foldchanges[global_indices])))
            rankings_gene_names.append(adata_comp.var_names[global_indices])
            rankings_gene_pvals.append(pvals[global_indices])
            rankings_gene_pvals_adj.append(pvals_adj[global_indices])
    elif method == 'wilcoxon':
        CONST_MAX_SIZE = 10000000
        means = np.zeros((n_groups, n_genes))
        vars = np.zeros((n_groups, n_genes))
        # initialize space for z-scores
        scores = np.zeros(n_genes)
        # First loop: Loop over all genes
        if reference != 'rest':
            for imask, mask in enumerate(groups_masks):
                means[imask], vars[imask] = get_mean_var(X[mask])  # for fold-change only
                if imask == ireference:
                    continue
                else:
                    mask_rest = groups_masks[ireference]
                ns_rest = np.where(mask_rest)[0].size
                mean_rest, var_rest = get_mean_var(X[mask_rest])  # for fold-change only

                if ns_rest <= 25 or ns[imask] <= 25:
                    log.info('Few observations in a group for '
                             'normal approximation (<=25). Lower test accuracy.')
                n_active = ns[imask]
                m_active = ns_rest

                # Now calculate gene expression ranking in chunkes:
                chunk = []
                # Calculate chunk frames
                n_genes_max_chunk = floor(CONST_MAX_SIZE / (n_active + m_active))
                if n_genes_max_chunk < n_genes:
                    chunk_index = n_genes_max_chunk
                    while chunk_index < n_genes:
                        chunk.append(chunk_index)
                        chunk_index = chunk_index + n_genes_max_chunk
                    chunk.append(n_genes)
                else:
                    chunk.append(n_genes)

                left = 0
                # Calculate rank sums for each chunk for the current mask
                for chunk_index, right in enumerate(chunk):
                    # Check if issparse is true: AnnData objects are currently sparse.csr or ndarray.
                    if issparse(X):
                        df1 = pd.DataFrame(data=X[mask, left:right].todense())
                        df2 = pd.DataFrame(data=X[mask_rest, left:right].todense(),
                                           index=np.arange(start=n_active, stop=n_active + m_active))
                    else:
                        df1 = pd.DataFrame(data=X[mask, left:right])
                        df2 = pd.DataFrame(data=X[mask_rest, left:right],
                                           index=np.arange(start=n_active, stop=n_active + m_active))
                    df1 = df1.append(df2)
                    ranks = df1.rank()
                    # sum up adjusted_ranks to calculate W_m,n
                    scores[left:right] = np.sum(ranks.loc[0:n_active, :])
                    left = right

                # https://en.wikipedia.org/wiki/Mann%E2%80%93Whitney_U_test
                # for large samples, U (score here) is approximately normally distributed
                scores = (scores - (n_active * (n_active + m_active + 1) / 2)) / sqrt(
                    (n_active * m_active * (n_active + m_active + 1) / 12))

                # for methylation, we are looking for reverse difference
                scores = scores * -1

                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores))

                if corr_method == 'benjamini-hochberg':
                    pvals[np.isnan(pvals)] = 1  # set Nan values to 1 to properly convert using Benhjamini Hochberg
                    _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
                elif corr_method == 'bonferroni':
                    pvals_adj = np.minimum(pvals * n_genes, 1.0)
                else:
                    raise ValueError

                # Fold change
                foldchanges = (means[imask] + 1e-9) / (mean_rest + 1e-9)  # add small value to remove 0's
                scores_sort = np.abs(scores) if rankby_abs else scores
                partition = np.argpartition(scores_sort, -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(scores_sort[partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_scores.append(scores[global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])
                rankings_gene_logfoldchanges.append(np.log2(foldchanges[global_indices]))
                rankings_gene_pvals.append(pvals[global_indices])
                rankings_gene_pvals_adj.append(pvals_adj[global_indices])

        # If no reference group exists, ranking needs only to be done once (full mask)
        else:
            scores = np.zeros((n_groups, n_genes))
            chunk = []
            n_cells = X.shape[0]
            n_genes_max_chunk = floor(CONST_MAX_SIZE / n_cells)
            if n_genes_max_chunk < n_genes:
                chunk_index = n_genes_max_chunk
                while chunk_index < n_genes:
                    chunk.append(chunk_index)
                    chunk_index = chunk_index + n_genes_max_chunk
                chunk.append(n_genes)
            else:
                chunk.append(n_genes)
            left = 0
            for chunk_index, right in enumerate(chunk):
                # Check if issparse is true
                if issparse(X):
                    df1 = pd.DataFrame(data=X[:, left:right].todense())
                else:
                    df1 = pd.DataFrame(data=X[:, left:right])
                ranks = df1.rank()
                # sum up adjusted_ranks to calculate W_m,n
                for imask, mask in enumerate(groups_masks):
                    scores[imask, left:right] = np.sum(ranks.loc[mask, :])
                left = right

            for imask, mask in enumerate(groups_masks):
                mask_rest = ~groups_masks[imask]
                means[imask], vars[imask] = get_mean_var(X[mask])  # for fold-change
                mean_rest, var_rest = get_mean_var(X[mask_rest])  # for fold-change

                scores[imask, :] = (scores[imask, :] - (ns[imask] * (n_cells + 1) / 2)) / sqrt(
                    (ns[imask] * (n_cells - ns[imask]) * (n_cells + 1) / 12))
                scores[np.isnan(scores)] = 0
                # for methylation, we are looking for reverse difference
                scores = scores * -1

                pvals = 2 * stats.distributions.norm.sf(np.abs(scores[imask, :]))

                if corr_method == 'benjamini-hochberg':
                    pvals[np.isnan(pvals)] = 1  # set Nan values to 1 to properly convert using Benhjamini Hochberg
                    _, pvals_adj, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')
                elif corr_method == 'bonferroni':
                    pvals_adj = np.minimum(pvals * n_genes, 1.0)
                else:
                    raise ValueError

                # Fold change
                foldchanges = (means[imask] + 1e-9) / (mean_rest + 1e-9)  # add small value to remove 0's
                scores_sort = np.abs(scores) if rankby_abs else scores
                partition = np.argpartition(scores_sort[imask, :], -n_genes_user)[-n_genes_user:]
                partial_indices = np.argsort(scores_sort[imask, partition])[::-1]
                global_indices = reference_indices[partition][partial_indices]
                rankings_gene_scores.append(scores[imask, global_indices])
                rankings_gene_names.append(adata_comp.var_names[global_indices])
                rankings_gene_logfoldchanges.append(np.log2(foldchanges[global_indices]))
                rankings_gene_pvals.append(pvals[global_indices])
                rankings_gene_pvals_adj.append(pvals_adj[global_indices])
    else:
        raise NotImplementedError(f'{method} is not implemented')

    groups_order_save = [str(g) for g in groups_order]
    if (reference != 'rest' and method != 'logreg') or (method == 'logreg' and len(groups) == 2):
        groups_order_save = [g for g in groups_order if g != reference]
    adata.uns[key_added]['scores'] = np.rec.fromarrays(
        [n for n in rankings_gene_scores],
        dtype=[(rn, 'float32') for rn in groups_order_save])
    adata.uns[key_added]['names'] = np.rec.fromarrays(
        [n for n in rankings_gene_names],
        dtype=[(rn, 'U50') for rn in groups_order_save])

    if method in {'t-test', 't-test_overestim_var', 'wilcoxon'}:
        adata.uns[key_added]['logfoldchanges'] = np.rec.fromarrays(
            [n for n in rankings_gene_logfoldchanges],
            dtype=[(rn, 'float32') for rn in groups_order_save])
        adata.uns[key_added]['pvals'] = np.rec.fromarrays(
            [n for n in rankings_gene_pvals],
            dtype=[(rn, 'float64') for rn in groups_order_save])
        adata.uns[key_added]['pvals_adj'] = np.rec.fromarrays(
            [n for n in rankings_gene_pvals_adj],
            dtype=[(rn, 'float64') for rn in groups_order_save])
    return adata if copy else None


def batch_correct_pc(adata, batch_series, correct=False,
                     n_components=30, sigma=25, alpha=0, knn=30,
                     **scanorama_kws):
    """
    Batch correction PCA based on scanorama

    Parameters
    ----------
    adata
        one major adata
    batch_series
        batch_series used for splitting adata
    correct
        if True, adata.X will be corrected inplace, otherwise only corrected PCs are added to adata.obsm['X_pca']
    dimred
        number of components in PCA
    sigma
        Correction smoothing parameter on Gaussian kernel.
    alpha
        Alignment score minimum cutoff.
    knn
        Number of nearest neighbors to use for matching.
    scanorama_kws
        Other Parameters passed to scanorama function
    Returns
    -------
    adata
    """
    try:
        import scanorama
    except ModuleNotFoundError as e:
        print('In order to use batch_correct_pc, you need to install scanorama, '
              'https://github.com/brianhie/scanorama')
        raise e

    scanorama_kws['dimred'] = n_components
    scanorama_kws['sigma'] = sigma
    scanorama_kws['alpha'] = alpha
    scanorama_kws['knn'] = knn

    adata.obs['batch'] = batch_series
    adata_list = []
    indexes = []
    for _, sub_df in adata.obs.groupby('batch'):
        adata_list.append(adata[sub_df.index, :])
        indexes.extend(sub_df.index.tolist())

    if correct:
        integrated, corrected = scanorama.correct_scanpy(adata_list, return_dimred=True,
                                                         **scanorama_kws)
        adata.X = np.vstack([ann.X.toarray() for ann in corrected])
    else:
        integrated = scanorama.integrate_scanpy(adata_list, **scanorama_kws)

    pca_df = pd.DataFrame(np.vstack(integrated), index=indexes).reindex(adata.obs_names)

    adata.obsm['X_pca'] = pca_df.values
    # TODO fill up other PC related parts same as sc.tl.pca
    return adata


def scanpy_umap(adata,
                min_dist=0.5,
                spread=1.0,
                n_components=2,
                maxiter=None,
                alpha=1.0,
                gamma=1.0,
                negative_sample_rate=5,
                init_pos='spectral',
                random_state=0,
                a=None,
                b=None):
    """Scanpy's umap function, take out here for getting umap but not adding to adata"""
    from scanpy.neighbors.umap.umap_ import find_ab_params, simplicial_set_embedding
    from scanpy.tools._utils import get_init_pos_from_paga
    if a is None or b is None:
        a, b = find_ab_params(spread, min_dist)
    else:
        a = a
        b = b
    if init_pos in adata.obsm.keys():
        init_coords = adata.obsm[init_pos]
    elif init_pos == 'paga':
        init_coords = get_init_pos_from_paga(adata, random_state=random_state)
    else:
        init_coords = init_pos
    from sklearn.utils import check_random_state
    random_state = check_random_state(random_state)
    n_epochs = maxiter
    _umap = simplicial_set_embedding(
        graph=adata.uns['neighbors']['connectivities'].tocoo(),
        n_components=n_components,
        initial_alpha=alpha,
        a=a,
        b=b,
        gamma=gamma,
        negative_sample_rate=negative_sample_rate,
        n_epochs=n_epochs,
        init=init_coords,
        random_state=random_state,
        verbose=False)
    return _umap


def run_embedding(adata, components=None, method='umap', method_kws=None):
    ALLOW_METHOD = {'umap', 'tsne'}

    if method_kws is None:
        method_kws = {}

    if components is None:
        components = [2]
    elif isinstance(components, int):
        components = [components]
    else:
        components = components

    method = method.lower()
    method_kws.pop('n_components', None)
    records = {}
    if method not in ALLOW_METHOD:
        raise ValueError(f'Method {method} not allowed in this function. Currently support {ALLOW_METHOD}.')
    else:
        _coord = None
        for component in components:
            if method == 'tsne':
                params_sklearn = dict(
                    perplexity=30,
                    random_state=0,
                    early_exaggeration=12,
                    learning_rate=1000,
                )
                params_sklearn.update(method_kws)

                from MulticoreTSNE import MulticoreTSNE as TSNE
                tsne = TSNE(n_components=component, **params_sklearn)
                # need to transform to float64 for MulticoreTSNE...
                _coord = tsne.fit_transform(adata.X.astype('float64'))
            elif method == 'umap':
                if component != 2:
                    method_kws['init_pos'] = 'spectral'
                _coord = scanpy_umap(adata, n_components=component, **method_kws)
            for dim in range(component):
                records[f'{method}_{component}_{dim}'] = _coord[:, dim]
        # put the last X_coord into adata
        if _coord is not None:
            adata.obsm[f'X_{method}'] = _coord
    total_df = pd.DataFrame(records, index=adata.obs_names)
    return total_df
