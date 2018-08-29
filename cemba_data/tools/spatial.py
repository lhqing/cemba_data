import h5py
import pandas as pd
import numpy as np
import holoviews as hv


def load_metainfo(hdf_path):
    # keys in allen_ish_grid.h5: ['dataset_meta', 'gene_meta', 'grid', 'grid_anno']
    with h5py.File(hdf_path) as allen_ish_grid:
        grid_anno = grid_3d(allen_ish_grid['grid_anno'].value)
        gene_meta = pd.DataFrame(allen_ish_grid['gene_meta'].value).set_index('id')
        gene_meta['gene_symbol'] = gene_meta['gene_symbol'].apply(lambda i: i.strip('*'))

        dataset_meta = pd.DataFrame(allen_ish_grid['dataset_meta'].value).set_index('section_data_set_id')
    return grid_anno, gene_meta, dataset_meta


def grid_3d(data):
    return data.reshape([58, 41, 67])


def grid_color_3d(data):
    return data.reshape([58, 41, 67, 3])


def left_brain_mask(data_3d):
    mask = np.zeros(data_3d.shape)
    sagittal_middle_index = int(data_3d.shape[0] / 2)
    mask[:sagittal_middle_index, :, :] = 1
    return mask


def get_roi_mask(anno_data, region_of_interest='all', descendant_map=None, left_half=False):
    if isinstance(region_of_interest, str) and region_of_interest == 'all':
        mask = anno_data != 0
    else:
        if isinstance(region_of_interest, str):
            region_of_interest = [region_of_interest]
        descendant_list = []
        for r in region_of_interest:
            descendant_list += descendant_map[r.lower()]
        mask = np.isin(anno_data, descendant_list)
    if left_half:  # only select left brain for sagittal dataset
        left_half_mask = left_brain_mask(anno_data)
        mask = mask * left_half_mask
    return mask.astype(bool)


def get_roi(anno_data, region_of_interest, descendant_map, left_half=False):
    mask = get_roi_mask(anno_data, region_of_interest, descendant_map, left_half)
    return anno_data * mask


def get_color(anno_data, color_map):
    anno_1d = np.ravel(anno_data)
    color_1d = np.array([color_map[i] for i in anno_1d])
    anno_color = color_1d.reshape([58, 41, 67, 3]) / 255
    return anno_color


def plot_anno(color_data, horizontal='max', sagittal='max', coronal='max'):
    non_white_pionts = np.argwhere(color_data.sum(axis=3) < 3)
    if sagittal == 'max':
        unique, counts = np.unique(non_white_pionts[:, 0], return_counts=True)
        sagittal = unique[counts.argmax()]
    if horizontal == 'max':
        unique, counts = np.unique(non_white_pionts[:, 1], return_counts=True)
        horizontal = unique[counts.argmax()]
    if coronal == 'max':
        unique, counts = np.unique(non_white_pionts[:, 2], return_counts=True)
        coronal = unique[counts.argmax()]

    sagittal_plot = hv.RGB(color_data[sagittal, :, :, :], bounds=True)
    horizontal_plot = hv.RGB(color_data[:, horizontal, :, :], bounds=True)
    coronal_plot = hv.RGB(np.swapaxes(color_data[:, :, coronal, :], 0, 1), bounds=True)
    holomap = hv.HoloMap({'sagittal': sagittal_plot,
                          'horizontal': horizontal_plot,
                          'coronal': coronal_plot}, kdims='plane')
    return holomap


def plot_grid_data(grid_data, horizontal=29, sagittal=20, coronal=33):
    if len(grid_data.shape) == 1:
        grid_data = grid_3d(grid_data.copy())

    sagittal_plot = hv.Image(grid_data[sagittal, :, :], bounds=True).options(cmap='viridis')
    horizontal_plot = hv.Image(grid_data[:, horizontal, :], bounds=True).options(cmap='viridis')
    coronal_plot = hv.Image(np.swapaxes(grid_data[:, :, coronal], 0, 1), bounds=True).options(cmap='viridis')
    holomap = hv.HoloMap({'sagittal': sagittal_plot,
                          'horizontal': horizontal_plot,
                          'coronal': coronal_plot}, kdims='plane')
    return holomap


def plot_gene(gene_grid, roi_color=None, horizontal=29, sagittal=20, coronal=33):
    gene_plots = plot_grid_data(gene_grid,
                                horizontal=horizontal,
                                sagittal=sagittal,
                                coronal=coronal)
    if roi_color is not None:
        anno_plots = plot_anno(roi_color,
                               horizontal=horizontal,
                               sagittal=sagittal,
                               coronal=coronal)
    else:
        anno_plots = None
    return gene_plots, anno_plots


def load_grid_data(hdf_path, col_select=None, as_df=False, columns=None):
    # keys in allen_ish_grid.h5: ['dataset_meta', 'gene_meta', 'grid', 'grid_anno']
    with h5py.File(hdf_path) as allen_ish_grid:
        if col_select is None:
            data = None
        elif col_select != 'all':
            data = allen_ish_grid['grid'][:, col_select]
        else:
            data = allen_ish_grid['grid'].value
    if as_df:
        data = pd.DataFrame(data, columns=columns)
    return data


def make_bins(length, bin_size=100):
    return [(i, min(i + bin_size, length))
            for i in range(0, length, bin_size)]


def dataset_qc(raw_data, valid_mask, above_threshold=5):
    _valid_mask = np.ravel(valid_mask)
    valid_data = raw_data[_valid_mask, :]
    mean = valid_data.mean(axis=0)
    above = (valid_data > above_threshold).sum(axis=0)
    return mean, above


def preprocess_data(raw_data, valid_mask):
    _data = raw_data[np.ravel(valid_mask), :]
    _data = np.log2(_data + 1)
    mean = _data.mean(axis=0)
    sd = _data.std(axis=0)
    _data = (_data - mean) / sd
    return _data


def add_df2group(h5_group, name, df, compression='gzip'):
    # save index into col
    df = df.copy()
    df.reset_index(drop=False, inplace=True)

    # judge dtype
    h5_dtype = []
    for col, col_dt in df.dtypes.iteritems():
        if col_dt == 'O':
            h5_dtype.append((col, h5py.special_dtype(vlen=str)))
        elif col_dt == 'int':
            h5_dtype.append((col, np.int32))
        elif col_dt == 'float':
            h5_dtype.append((col, np.float32))
        else:
            h5_dtype.append((col, col_dt))
    # h5 compound dtype
    h5_dtype = np.dtype(h5_dtype)
    h5_data = [tuple(row) for i, row in df.iterrows()]

    h5_group.require_dataset(name, shape=(df.shape[0],), dtype=h5_dtype, compression=compression)
    h5_group[name][...] = h5_data
    return
