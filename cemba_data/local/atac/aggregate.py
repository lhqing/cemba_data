import anndata
import pandas as pd
import scipy.sparse as ss


def aggregate_adata(cell_anno_df,
                    sample_adata_path_dict, sample_key_col,
                    cell_type_select, cell_type_col='MajorType'):
    cell_type_anno = cell_anno_df[cell_anno_df[cell_type_col].isin(cell_type_select)]
    print(f'Total cell: {cell_type_anno.shape[0]}')
    adata_x_list = []
    adata_obs_list = []
    adata = None
    for sample, path in sample_adata_path_dict.items():
        print(f'Opening {sample}')
        _type_cell_anno = cell_type_anno[cell_type_anno[sample_key_col] == sample]
        adata = anndata.read_h5ad(path)
        adata = adata[_type_cell_anno.index, :]
        for col_name, col in _type_cell_anno.iteritems():
            adata.obs[col_name] = col
        adata_x_list.append(adata.X)
        adata_obs_list.append(adata.obs)
    total_x = ss.vstack(adata_x_list)
    total_cell_obs = pd.concat(adata_obs_list, sort=True)
    total_cell_obs.index = total_cell_obs['Sample'] + '-' + total_cell_obs.index

    total_adata = anndata.AnnData(X=total_x.tocsr(),
                                  obs=total_cell_obs,
                                  var=adata.var)
    return total_adata


def aggregate_snap():
    """Read snap to adata, then select the cells, aggregate and return adata"""
    raise NotImplementedError
