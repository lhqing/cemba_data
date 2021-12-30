import pandas as pd
import h5py
import numpy as np
import pyBigWig
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from ALLCools.utilities import parse_chrom_size


def dump_single_snap(output_prefix, snap_path, barcode_to_cluster):
    cluster_handle_dict = {cluster: open(f'{output_prefix}_{cluster}.bed', 'w')
                           for cluster in barcode_to_cluster.unique()}
    cluster_n_frags = defaultdict(int)
    with h5py.File(snap_path, 'r') as snap:
        barcode = pd.Series(np.array(snap['/BD/name']).astype(str))
        barcode_length = np.array(snap['/FM/barcodeLen'])
        barcode_pos = np.array(snap['/FM/barcodePos']) - 1
        barcode_table = pd.DataFrame({
            'barcode': barcode,
            'length': barcode_length,
            'start': barcode_pos
        })
        barcode_table['end'] = barcode_table['start'] + barcode_table['length']
        barcode_table['cluster'] = barcode.map(barcode_to_cluster)
        barcode_table.dropna(inplace=True)

        for i, (_, row) in enumerate(barcode_table.iterrows()):
            *_, start, end, cluster = row
            frag_chrom = pd.Series(
                np.array(snap['/FM/fragChrom'][start:end]).astype(str))
            frag_start = pd.Series(np.array(snap['/FM/fragStart'][start:end]))
            frag_length = pd.Series(np.array(snap['/FM/fragLen'][start:end]))
            frag_table = pd.DataFrame({
                'chrom': frag_chrom,
                'start': frag_start,
                'end': frag_start + frag_length
            })
            cluster_handle_dict[cluster].write(frag_table.to_csv(sep='\t',
                                                                 header=None,
                                                                 index=None))
            cluster_n_frags[cluster] += frag_table.shape[0]

    for h in cluster_handle_dict.values():
        h.close()
    return cluster_n_frags


def fragments_to_bigwig(output_prefix, cluster, chrom_size_path, bw_bin_size=10, scale=None):
    # concat bed
    subprocess.run(
        f'cat {output_prefix}_*_{cluster}.bed > {output_prefix}_{cluster}.bed && '
        f'rm -f {output_prefix}_*_{cluster}.bed',
        shell=True,
        check=True)
    # sort bed
    subprocess.run(f'sort -k 1,1 -k2,2n {output_prefix}_{cluster}.bed '
                   f'> {output_prefix}_{cluster}.sorted.bed && '
                   f'rm -f {output_prefix}_{cluster}.bed',
                   shell=True,
                   check=True)
    # bed to bedgraph
    subprocess.run(f'bedtools genomecov -i {output_prefix}_{cluster}.sorted.bed '
                   f'-g {chrom_size_path} -bg > {output_prefix}_{cluster}.bedgraph',
                   shell=True,
                   check=True)
    bg_path = f'{output_prefix}_{cluster}.bedgraph'
    bw_path = f'{output_prefix}_{cluster}.bw'

    bg_iter = pd.read_csv(bg_path,
                          sep='\t',
                          header=None,
                          names=['chrom', 'start', 'end', 'count'],
                          chunksize=100000)
    total_wigs = []
    for bg in bg_iter:
        bg['start'] = bg['start'] // bw_bin_size
        wig_values = bg.groupby(['chrom',
                                 'start'])['count'].mean().reset_index()
        total_wigs.append(wig_values)
    total_wigs = pd.concat(total_wigs)
    total_wigs = total_wigs.groupby(['chrom',
                                     'start'])['count'].mean().reset_index()

    chrom_sizes = parse_chrom_size(chrom_size_path)
    chrom_sizes_list = [(k, v) for k, v in chrom_sizes.items()]

    with pyBigWig.open(bw_path, 'w') as bw_out:
        bw_out.addHeader(chrom_sizes_list)
        for chrom in chrom_sizes.keys():
            chrom_df = total_wigs[total_wigs['chrom'] == chrom].sort_values(
                'start')
            if chrom_df.shape[0] == 0:
                continue
            if scale is None:
                values = chrom_df['count'].astype(float)
            else:
                values = chrom_df['count'].astype(float) / scale
            bw_out.addEntries(chrom, (chrom_df['start'] * bw_bin_size).tolist(),
                              values=values.tolist(),
                              span=bw_bin_size)

    # remove bed graph
    subprocess.run(f'rm -f {bg_path}', shell=True, check=True)
    # gzip bed
    subprocess.run(f'gzip {output_prefix}_{cluster}.sorted.bed', shell=True, check=True)
    return


def atac_bulk(cluster_table, snap_table, output_prefix, bw_bin_size=10, cpu=1):
    cluster_table = pd.read_csv(cluster_table, index_col=None, header=None,
                                names=['sample', 'barcode', 'cluster'])
    snap_table = pd.read_csv(snap_table, index_col=0, header=None, squeeze=True)

    sample_clusters = {s: c.set_index('barcode')['cluster'].to_dict()
                       for s, c in cluster_table.groupby('sample')}
    for s in sample_clusters.keys():
        c = sample_clusters[s]
        c.index = c.index.str.split('-').str[1]

    # dump fragments
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for sample, clusters in sample_clusters.items():
            _output_prefix = f'{output_prefix}_{sample}'
            f = exe.submit(dump_single_snap,
                           output_prefix=_output_prefix,
                           snap_path=snap_table[sample],
                           barcode_to_cluster=clusters)
            futures[f] = sample

        total_n_frags = defaultdict(int)
        for f in as_completed(futures):
            sample = futures[f]
            print(f'{sample} returned')
            cluster_n_frags = f.result()
            for k, v in cluster_n_frags.items():
                total_n_frags[k] += v

    million_reads_scaler = {k: v / 1000000 for k, v in total_n_frags.items()}
    total_clusters = pd.concat(sample_clusters.values()).unique()
    # generate bed and bigwig
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for cluster in total_clusters:
            f = exe.submit(fragments_to_bigwig,
                           output_prefix,
                           cluster,
                           bw_bin_size=bw_bin_size,
                           scale=million_reads_scaler['ASC'])
            futures[f] = cluster

        for f in as_completed(futures):
            cluster = futures[f]
            print(f'{cluster} returned')
            f.result()
