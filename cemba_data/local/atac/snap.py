import pandas as pd
import numpy as np
import scipy.sparse as ss
from anndata import AnnData
import pybedtools
import h5py


def read_snap(file_path, bin_size=5000):
    with h5py.File(file_path) as f:
        data = f[f'/AM/{bin_size}/count'].value
        idx = f[f'/AM/{bin_size}/idx'].value - 1  # R is 1 based, python is 0 based
        idy = f[f'/AM/{bin_size}/idy'].value - 1  # R is 1 based, python is 0 based

        bin_chrom = f[f'/AM/{bin_size}/binChrom'].value.astype(str)
        bin_start = f[f'/AM/{bin_size}/binStart'].value - 1  # 0 based bed format
        bin_end = f[f'/AM/{bin_size}/binStart'].value - 1 + bin_size  # 0 based bed format

        bin_id = np.core.defchararray.add(np.core.defchararray.add(bin_chrom, '-'),
                                          bin_start.astype(str))

        cell_barcode = f['/BD/name'].value.astype(str)
        cell_id = [f'ATAC_{i}' for i in range(cell_barcode.size)]
        cell_meta = pd.DataFrame([], index=cell_id)
        cell_meta['barcode'] = cell_barcode
        for name in f['/BD/'].keys():
            if name == 'name':
                continue
            cell_meta[name] = f[f'/BD/{name}'].value
        data = ss.coo_matrix((data, (idx, idy)), shape=(cell_barcode.size, bin_id.size)).tocsc()

    adata = AnnData(X=data,
                    obs=cell_meta,
                    var=pd.DataFrame({'chrom': bin_chrom,
                                      'start': bin_start,
                                      'end': bin_end},
                                     index=pd.Index(bin_id, name='chrom5k')))
    return adata


def reshape_matrix(adata, bed_path, chrom_size_path, slop=5000, force_run=False):
    """
    Reshape adata.X columns into any region list, based on bedtools intersect.
    Return a CSC matrix. for 60K genes and 60K cells, take ~90s.
    Note: This function internally generate a full matrix, don't use this for small regions.

    Parameters
    ----------
    adata
    bed_path
    chrom_size_path
    slop
    force_run

    Returns
    -------

    """
    # TODO return adata
    region_bed = pd.read_table(bed_path, index_col=-1,
                               header=None, names=['chrom', 'start', 'end', 'id'])
    if (not force_run) and (region_bed.shape[0] * adata.shape(0) > 1e10):
        raise ValueError('This function internally generate a full matrix, '
                         'do not use this for too many regions.')
    else:
        print('Note: This function internally generate a full matrix, '
              'do not use this for too many regions.')
    region_bt = pybedtools.BedTool.from_dataframe(region_bed.reset_index().iloc[:, [1, 2, 3, 0]])
    # for bin_id, instead of using text, use number
    bin_bt = pybedtools.BedTool.from_dataframe(adata.var
                                               .reset_index()
                                               .reset_index()
                                               .iloc[:, [2, 3, 4, 0, 1]])
    if slop is None:
        intersected_bin = region_bt.intersect(bin_bt, wa=True, wb=True)
    else:
        intersected_bin = region_bt.slop(b=slop, g=chrom_size_path) \
            .intersect(bin_bt, wa=True, wb=True)
    result = intersected_bin.to_dataframe()
    result.columns = ['chrom', 'start', 'end', 'region_id',
                      'bin_chrom', 'bin_start', 'bin_end', 'bin_num', 'bin_name']

    region_data = np.zeros((adata.obs_names.size,
                            result['region_id'].unique().size),
                           dtype=np.uint16)
    csc_data = adata.X.tocsc()
    for col_idx, (region, sub_df) in enumerate(result.sort_values('bin_num').groupby('region_id')):
        cell_region_count = csc_data[:, slice(sub_df.iat[0, 7],
                                              sub_df.iat[-1, 7] + 1)].sum(axis=1)
        cell_region_count = cell_region_count.A1
        region_data[:, col_idx] = cell_region_count
    sparse_result = ss.csc_matrix(region_data)
    return sparse_result


def reshape_matrix_fix_step(adata, window, step):
    """
    Reshape adata.X columns based on fix window size and steps.
    Window and step refer to matrix columns, not any genome region size.
    """
    batch = 3000

    chrom, start = zip(*adata.var_names.str.split('-'))
    cur_chrom = None
    chrom_start_cols = []
    unique_chrom_list = []
    for i, c in enumerate(chrom):
        if c != cur_chrom:
            chrom_start_cols.append(i)
            unique_chrom_list.append(c)
        cur_chrom = c
    chrom_start_cols.append(len(chrom))
    chrom_slices = []
    for i, c in enumerate(unique_chrom_list):
        start, end = chrom_start_cols[i:i + 2]
        if end - start > window:
            chrom_slices.append(slice(start, end))

    rows, _ = adata.X.shape
    csr_data = adata.X.tocsr()
    chunk_generator = (csr_data[i:i + batch, :] for i in range(0, rows, batch))

    csr_results = []
    n = 0
    for chunk in chunk_generator:
        sub_data = chunk.toarray()
        n += 1
        print(n)
        chrom_results = []
        for chrom_slice in chrom_slices:
            chrom_data = sub_data[:, chrom_slice]
            data_list = [chrom_data[:, slice(w, chrom_data.shape[1] - window + w, step)]
                         for w in range(window)]
            sum_of_data = np.sum(data_list, axis=0)
            chrom_results.append(sum_of_data)
        chunk_total = np.hstack(chrom_results)
        csr_results.append(ss.csr_matrix(chunk_total))
    total_result = ss.vstack(csr_results, dtype=np.uint32)

    # slice the index with same strategy
    total_chrom_index = []
    for chrom_slice in chrom_slices:
        col_index = adata.var_names[chrom_slice]
        chrom_col_index = col_index[slice(0, col_index.size - window, step)].tolist()
        total_chrom_index += chrom_col_index
    total_chrom_index = pd.Index(total_chrom_index)

    result_adata = AnnData(X=total_result,
                           obs=adata.obs,
                           var=pd.DataFrame([], index=total_chrom_index))

    return result_adata


def split_bam():
    # TODO: given a cluster assignment, split bam file
    return


def merge_cell():
    # TODO: generate meta cell based on a given embedding or cell index
    return
