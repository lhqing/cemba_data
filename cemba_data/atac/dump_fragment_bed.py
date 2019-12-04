import collections
import pathlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import lru_cache

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as ss
from ALLCools._open import open_gz
from anndata import AnnData


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


def read_snap_gene(file_path):
    with h5py.File(file_path) as f:
        data = f[f'/GM/count'].value
        idx = f[f'/GM/idx'].value - 1  # R is 1 based, python is 0 based
        idy = f[f'/GM/idy'].value - 1  # R is 1 based, python is 0 based

        gene_id = f[f'/GM/name'].value.astype(str)

        cell_barcode = f['/BD/name'].value.astype(str)
        cell_id = [f'ATAC_{i}' for i in range(cell_barcode.size)]
        cell_meta = pd.DataFrame([], index=cell_id)
        cell_meta['barcode'] = cell_barcode
        for name in f['/BD/'].keys():
            if name == 'name':
                continue
            cell_meta[name] = f[f'/BD/{name}'].value
        data = ss.coo_matrix((data, (idx, idy)), shape=(cell_barcode.size, gene_id.size)).tocsc()

    adata = AnnData(X=data,
                    obs=cell_meta,
                    var=pd.DataFrame([],
                                     index=pd.Index(gene_id, name='gene')))
    return adata


@lru_cache(maxsize=100)
def _get_barcode_dict(snap_path):
    with h5py.File(snap_path, 'r') as f:
        barcode_dict = collections.OrderedDict()
        i = 0
        for item in f["BD/name"]:
            item = item.decode()
            barcode_dict[item] = i
            i = i + 1
    return barcode_dict


def _get_frags_iter(f, barcodes, barcode_dict, barcode_chunk=50):
    frag_list = []
    count = 0
    for barcode in barcodes:
        try:
            barcode_id = barcode_dict[barcode]
        except KeyError:
            print(f'barcode {barcode} not found in snap file')
            continue
        barcode_pos = f["FM"]["barcodePos"][barcode_id] - 1
        barcode_len = f["FM"]["barcodeLen"][barcode_id]
        _chroms = [item.decode() for item in f["FM"]["fragChrom"][barcode_pos:(barcode_pos + barcode_len)]]
        _start = f["FM"]["fragStart"][barcode_pos:(barcode_pos + barcode_len)]
        _len = f["FM"]["fragLen"][barcode_pos:(barcode_pos + barcode_len)]
        frag_list = frag_list + list(zip(_chroms, _start, _start + _len, [1] * len(_chroms)))

        count += 1
        if count % barcode_chunk == 0:
            yield frag_list
            frag_list = []
    yield frag_list


def _dump_frags_single_snap(snap_path, out_f, barcodes):
    total_frags = 0
    barcode_dict = _get_barcode_dict(snap_path)

    with h5py.File(snap_path, 'r') as f:
        for frag_chunk in _get_frags_iter(f, barcodes, barcode_dict, barcode_chunk=50):
            total_frags += len(frag_chunk)
            for frag in frag_chunk:
                out_f.write("\t".join(map(str, frag)) + "\n")
    return total_frags


def _dump_frags_single_cluster(output_path, cluster_df, sample_snap_dict):
    with open_gz(output_path, 'w') as f:
        # cluster_df col 0 is sample, col 1 is barcode
        for sample, cluster_sample_df in cluster_df.groupby(cluster_df.iloc[:, 0]):
            snap_path = sample_snap_dict[sample]
            barcodes = cluster_sample_df.iloc[:, 1].tolist()
            _dump_frags_single_snap(snap_path=snap_path, out_f=f, barcodes=barcodes)
    return


def dump_frags(cell_group_path, output_dir_path, sample_snap_path, cpu=1):
    """
    Extract fragment bed file for levels of clustering results from snap file.
    """
    sample_snap_dict = pd.read_csv(sample_snap_path, index_col=0, header=None, squeeze=True).to_dict()
    cell_group_table = pd.read_csv(cell_group_path, index_col=None)
    output_dir = pathlib.Path(output_dir_path).absolute()
    future_list = []
    output_path_list = []
    with ProcessPoolExecutor(cpu) as executor:
        for col in cell_group_table.columns[2:]:
            col_dir = output_dir / col
            col_dir.mkdir(exist_ok=True)

            for cluster, cluster_df in cell_group_table.groupby(col):
                output_path = col_dir / f'{cluster}.bed.gz'
                future = executor.submit(_dump_frags_single_cluster,
                                         output_path=output_path,
                                         cluster_df=cluster_df,
                                         sample_snap_dict=sample_snap_dict)
                future_list.append(future)
                output_path_list.append(output_path)
    for future in as_completed(future_list):
        future.result()
    return output_path_list
