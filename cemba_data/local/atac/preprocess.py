import numpy as np


def filter_atac_cell(adata,
                     fragment_num_upper=2000,
                     mito_ratio_lower=0.1,
                     umap_ratio_upper=0.8,
                     dup_ratio_lower=0.6,
                     pair_ratio_upper=0.95):
    adata.obs['fragment_num'] = adata.obs['TN']
    adata.obs['mito_ratio'] = adata.obs['CM'] / adata.obs['UQ']
    adata.obs['umap_ratio'] = adata.obs['UM'] / adata.obs['TN']
    adata.obs['dup_ratio'] = 1 - adata.obs['UQ'] / adata.obs['PP']
    adata.obs['pair_ratio'] = adata.obs['PP'] / adata.obs['UM']

    filtered_obs = adata.obs.query(f'(fragment_num > {fragment_num_upper}) & '
                                   f'(mito_ratio < {mito_ratio_lower}) & '
                                   f'(umap_ratio > {umap_ratio_upper}) & '
                                   f'(dup_ratio < {dup_ratio_lower}) & '
                                   f'(pair_ratio > {pair_ratio_upper})')
    adata = adata[filtered_obs.index, :]
    return adata


def filter_bins(adata, zscore_cutoff=2):
    x = adata.X
    bin_sum = x.sum(axis=0)
    lg1p_bin_sum = np.log1p(bin_sum)
    _mean = lg1p_bin_sum[lg1p_bin_sum > 0].mean()
    _sd = lg1p_bin_sum[lg1p_bin_sum > 0].std()
    norm_lg1p_bin_sum = (lg1p_bin_sum - _mean) / _sd
    bool_index = (norm_lg1p_bin_sum > -abs(zscore_cutoff)) & (norm_lg1p_bin_sum < abs(zscore_cutoff))
    _adata = adata[:, np.array(bool_index).reshape(bool_index.size)]
    return _adata


def scale_jaccard(jaccard_matrix,
                  with_row_mean=True, with_row_std=True, clip_value=5,
                  with_col_mean=True, with_col_std=False):
    clip_value = abs(clip_value)

    # scale row
    from sklearn.preprocessing import StandardScaler
    scaler = StandardScaler(with_mean=with_row_mean, with_std=with_row_std)
    jaccard_matrix = scaler.fit_transform(jaccard_matrix.T).T
    jaccard_matrix[jaccard_matrix > clip_value] = clip_value
    jaccard_matrix[jaccard_matrix < -clip_value] = -clip_value

    # scale col
    scaler = StandardScaler(with_mean=with_col_mean, with_std=with_col_std)
    jaccard_matrix = scaler.fit_transform(jaccard_matrix)
    return jaccard_matrix


def pca(adata, n_components=25, whiten=True, random_state=0, pca_kws=None):
    from sklearn.decomposition import PCA
    if pca_kws is None:
        pca_kws = {}
    pca = PCA(n_components=n_components, whiten=whiten,
              random_state=random_state, **pca_kws)
    pc = pca.fit_transform(adata.X)

    # weight pc by SGV
    # jaccard_pc = jaccard_pc @ np.diag(np.sqrt(pca.singular_values_))

    adata.obsm['X_pca'] = pc
    adata.varm['PCs'] = pca.components_.T
    adata.uns['pca'] = {}
    adata.uns['pca']['variance'] = pca.explained_variance_
    adata.uns['pca']['variance_ratio'] = pca.explained_variance_ratio_
    return adata
