import pandas as pd
import numpy as np
import scipy.sparse as ss
from anndata import AnnData
import pybedtools
import h5py
import subprocess
import pathlib


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
    Note: This function internally generate a full matrix, don't use this for >100k regions.

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
    region_bed = pd.read_table(bed_path, index_col=-1,
                               header=None, names=['chrom', 'start', 'end', 'id'])
    if (not force_run) and (region_bed.shape[0] * adata.shape[0] > 1e10):
        raise ValueError('This function internally generate a full matrix, '
                         'do not use this for too many regions.')
    else:
        print('Note: This function internally generate a full matrix, '
              'do not use this for too many regions.')
    region_bt = pybedtools.BedTool.from_dataframe(region_bed.reset_index().iloc[:, [1, 2, 3, 0]])
    # for bin_id, instead of using text, use number
    bin_bt = pybedtools.BedTool.from_dataframe(adata.var
                                               .reset_index()  # region_id
                                               .reset_index()  # int id
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
    regions = []
    for col_idx, (region, sub_df) in enumerate(result.sort_values('bin_num').groupby('region_id')):
        regions.append(region)
        cell_region_count = csc_data[:, slice(sub_df.iat[0, 7],
                                              sub_df.iat[-1, 7] + 1)].sum(axis=1).A1
        region_data[:, col_idx] = cell_region_count
    sparse_result = ss.csc_matrix(region_data)

    result_adata = AnnData(X=sparse_result,
                           obs=adata.obs.copy(),
                           var=pd.DataFrame([], index=pd.Index(regions)))
    return result_adata


# this is wrong, SNAP matrix feature is not full, bins without any signal in all cell are dropped.
"""
def reshape_matrix_fix_step(adata, window, step):
    \"""
    Reshape adata.X columns based on fix window size and steps.
    Window and step refer to matrix columns, not any genome region size.
    \"""
    batch = 3000

    # prepare chrom_slices to separate chromosome columns
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
        chrom_results = []
        for chrom_slice in chrom_slices:
            chrom_data = sub_data[:, chrom_slice]
            # slice the data based on window and step, sum to get result
            # this slice dropped the last several bins in each chromosome, which is usually OK
            # TODO: deal with last several bins to make this perfect, remember to change index slice too
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
                           obs=adata.obs.copy(),
                           var=pd.DataFrame([], index=total_chrom_index))
    return result_adata
"""


def split_bam(bam_path, out_prefix, cell_to_cluster, out_dir=None,
              keep_unknown=False, cpu=10, mapq_cutoff=10):
    """

    Parameters
    ----------
    bam_path

    out_prefix

    cell_to_cluster
        A dict where the key is cell barcode, value is cluster assignment
    out_dir

    keep_unknown

    cpu

    mapq_cutoff

    Returns
    -------

    """
    bam_path = pathlib.Path(bam_path).absolute()
    # get bam header text
    headers = subprocess.run(['samtools', 'view', '-H', str(bam_path)],
                             stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf8').stdout

    # set up out dir
    if out_dir is None:
        out_dir = bam_path.parent
    else:
        out_dir = pathlib.Path(out_dir).absolute()
    # set up output handle for each cluster
    unique_clusters = set(cell_to_cluster.values())
    if keep_unknown:
        unique_clusters.add('unknown')
    cluster_handle_dict = {}
    for cluster in unique_clusters:
        out_handle = open(out_dir / f'{out_prefix}.{cluster}.sam', 'w')
        out_handle.write(headers)
        cluster_handle_dict[cluster] = out_handle

    # read and split bam file
    read_popen = subprocess.Popen(['samtools', 'view', str(bam_path)],
                                  stdout=subprocess.PIPE, encoding='utf8')
    for line in read_popen.stdout:
        cell_id = line.split(':')[0]
        try:
            cluster_handle_dict[cell_to_cluster[cell_id]].write(line)
        except KeyError:
            if keep_unknown:
                cluster_handle_dict['unknown'].write(line)
            else:
                continue
    for k, v in cluster_handle_dict.items():
        v.close()

    # transfer all split sam into bam
    for cluster in unique_clusters:
        out_sam = str(out_dir / f'{out_prefix}.{cluster}.sam')
        out_bam = str(out_dir / f'{out_prefix}.{cluster}.bam')
        sort_bam = str(out_dir / f'{out_prefix}.{cluster}.sort.bam')
        dedup_bam = str(out_dir / f'{out_prefix}.{cluster}.dedup.bam')
        dedup_matrix = dedup_bam + '.matrix'
        cmd = f'samtools view -q {mapq_cutoff} -f 2 {out_sam} -o {out_bam}'
        subprocess.run(cmd.split(' '), check=True)
        cmd = f'samtools sort -@ {cpu} {out_bam} -o {sort_bam}'
        subprocess.run(cmd.split(' '), check=True)
        subprocess.run(['rm', '-f', out_sam])
        dedup_cmd = f'picard MarkDuplicates I={sort_bam} O={dedup_bam} ' \
                    f'M={dedup_matrix} REMOVE_DUPLICATES=true'
        subprocess.run(dedup_cmd.split(' '), check=True)
        cmd = ['samtools', 'index', str(dedup_bam)]
        subprocess.run(cmd, stderr=subprocess.PIPE)
        cmd = ['bamCoverage', '-b', str(dedup_bam), '-o',
               str(dedup_bam)[:-3] + 'bigwig', '--binSize', '50',
               '-p', '1', '--normalizeUsing', 'CPM']
        subprocess.run(cmd, stderr=subprocess.PIPE)
    return


def merge_cell(adata, embedding_data, target=5000, n_neighbors=30, max_dist=0.3, chunksize=500, seed=0):
    from sklearn.neighbors import LocalOutlierFactor

    # Calculate neighbor and density
    clf = LocalOutlierFactor(n_neighbors=n_neighbors, algorithm='auto',
                             leaf_size=30, metric='minkowski',
                             p=2, metric_params=None, contamination='auto', novelty=False, n_jobs=30)
    clf.fit(embedding_data)
    scores = clf.negative_outlier_factor_
    scores = (scores - scores.min()) / (scores.max() - scores.min())
    probability = scores / scores.sum()
    # density weighted cell selection
    np.random.seed(seed)
    cell_selection = np.random.choice(range(adata.obs_names.size), size=target,
                                      replace=False, p=probability)
    original_name = adata.obs_names[cell_selection]

    # get neighbors
    # TODO: try more fancy neighbor select
    distance, points = clf.kneighbors(embedding_data)
    distance = distance[cell_selection, :]
    points = points[cell_selection, :]
    distance_mask = distance < max_dist

    # calculate neighbor sum
    csr_data = adata.X.tocsr()
    n = 0
    chunk_result = np.zeros((chunksize, adata.shape[1]), dtype=np.uint32)
    csr_chunk_results = []
    for _dist, _point in zip(distance_mask, points):
        _use_point = _point[_dist]
        point_neigh_sum = csr_data[_use_point, :].sum(axis=0, dtype=np.uint32).A1
        chunk_result[n, :] = point_neigh_sum
        n += 1
        if n % chunksize == 0:
            n = 0
            csr_chunk_results.append(ss.csr_matrix(chunk_result))
            chunk_result = np.zeros((chunksize, adata.shape[1]), dtype=np.uint32)
    csr_result = ss.vstack(csr_chunk_results)

    cell_id = [f'METACELL_{i}' for i in range(csr_result.shape[0])]
    meta_adata = AnnData(X=csr_result,
                         obs=pd.DataFrame({'center_raw_id': original_name},
                                          index=pd.Index(cell_id)),
                         var=adata.var.copy())
    return meta_adata
