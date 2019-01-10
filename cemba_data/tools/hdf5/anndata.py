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


def highly_variable_feature(
        adata,
        min_disp=None, max_disp=None,
        min_mean=None, max_mean=None,
        n_top_feature=None,
        n_bins=50,
        no_info_cutoff=0.02):
    """
    Adapted from Scanpy, see license above
    """
    # TODO: add cov in bin average, basicaly a 3D scale,
    # RNA is only scaled by mean, but methylation is scaled by both mean and cov
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
    mean, var = get_mean_var(total_X)
    # calculate total mean var 1st, filter out those region that have no coverage. (mean -> 1, var -> 0)
    use_var = adata.var.index[(mean - 1) ** 2 + var > no_info_cutoff]
    X = adata[:, use_var].X
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
    other_mean, other_var = get_mean_var(other_var)
    other_df = pd.DataFrame({
        'mean': other_mean,
        'dispersion': df['dispersion'].min(),
        'dispersion_norm': df['dispersion_norm'].min(),
        'gene_subset': False,
    }, index=adata.var.index[~adata.var.index.isin(use_var)])
    df = pd.concat([df, other_df], sort=True).reindex(adata.var.index)
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
        data_type='mc',
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
    # TODO add other tests, and check out https://github.com/theislab/diffxpy and other info mentioned here:
    # https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html (marker finding part)

    log.info('ranking genes')
    avail_methods = {'t-test', 't-test_overestim_var'}
    if method not in avail_methods:
        raise ValueError('Method must be one of {}.'.format(avail_methods))

    avail_corr = {'benjamini-hochberg', 'bonferroni'}
    if corr_method not in avail_corr:
        raise ValueError('Correction method must be one of {}.'.format(avail_corr))

    adata = adata.copy() if copy else adata
    # sanitize anndata
    # adata._sanitize()

    # for clarity, rename variable
    groups_order = groups
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

            if data_type == 'mc':
                # for methylation, we reverse the score, this doesn't impact p-value, but impact score sorting.
                # because we are looking for hypo-methylation
                scores = (mean_rest - means[igroup]) / denominator
            else:
                scores = (means[igroup] - mean_rest) / denominator

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
    else:
        raise NotImplementedError(f'{method} is not implemented, buy me a coffee I can implement')

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
    for _, sub_df in adata.obs.groupby('batch'):
        adata_list.append(adata[sub_df.index, :])

    if correct:
        integrated, corrected = scanorama.correct_scanpy(adata_list, **scanorama_kws)
        adata.X = np.vstack([ann.X.toarray() for ann in corrected])
    else:
        integrated = scanorama.integrate_scanpy(adata_list, **scanorama_kws)
    adata.obsm['X_pca'] = np.vstack(integrated)

    # TODO fill up other PC related parts same as sc.tl.pca
    return adata
