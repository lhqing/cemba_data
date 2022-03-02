import pysam
import pandas as pd
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import os


def merge_single_bam(bam_path, cell_id_to_cluster, output_prefix, header_dict):
    header = pysam.AlignmentHeader.from_dict(header_dict)
    clusters = set(cell_id_to_cluster.values())
    cluster_read_counts = {c: 0 for c in clusters}

    # write reads by cluster
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        # open BAM handles for each cluster
        cluster_handles = {}
        for cluster in clusters:
            cluster_handles[cluster] = pysam.AlignmentFile(
                f'{output_prefix}_{cluster}.bam', "wb", header=header)

        for read in bam_file:
            cell_id = read.get_tag('RG')
            try:
                cluster = cell_id_to_cluster[cell_id]
                # this removes RG tag
                read.set_tag('RG', None)
                cluster_handles[cluster].write(read)
                cluster_read_counts[cluster] += 1
            except KeyError:
                continue

    # close handles
    for handle in cluster_handles.values():
        handle.close()

    # delete empty out_bam
    for cluster, count in cluster_read_counts.items():
        bam_path = f'{output_prefix}_{cluster}.bam'
        if count == 0:
            subprocess.run(['rm', '-rf', bam_path])
    return cluster_read_counts


def merge_mct_cluster_bam(cell_id_to_cluster_path,
                          bam_list_path,
                          output_prefix,
                          cpu=10):
    cell_id_to_cluster = pd.read_csv(
        cell_id_to_cluster_path,
        index_col=0,
        header=None,
        squeeze=True).to_dict()
    bam_paths = pd.read_csv(bam_list_path, header=None, squeeze=True).tolist()

    # get header
    with pysam.AlignmentFile(bam_paths[0]) as bam:
        header_dict = bam.header.as_dict()
        # remove cell specific info
        keys_to_delete = ['PG', 'RG', 'CO']
        for k in keys_to_delete:
            if k in header_dict:
                del header_dict[k]

    clusters = set(cell_id_to_cluster.values())
    total_cluster_read_counts = {c: 0 for c in clusters}

    # merge single bam files
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for i, path in enumerate(bam_paths):
            f = exe.submit(merge_single_bam,
                           bam_path=path,
                           cell_id_to_cluster=cell_id_to_cluster,
                           output_prefix=f'{output_prefix}{i:06d}',
                           header_dict=header_dict)
            futures[f] = path

        for f in as_completed(futures):
            cluster_read_counts = f.result()
            for k, v in cluster_read_counts.items():
                total_cluster_read_counts[k] += v

    # merge cluster bam files
    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for cluster in clusters:
            chunk_paths = list(glob.glob(f'{output_prefix}*_{cluster}.bam'))
            if len(chunk_paths) == 0:
                continue
            merge_cmd = f'samtools merge --no-PG -c -o {output_prefix}_{cluster}.bam ' \
                        f'{output_prefix}*_{cluster}.bam && ' \
                        f'samtools index {output_prefix}_{cluster}.bam'
            f = exe.submit(subprocess.run,
                           merge_cmd,
                           shell=True,
                           check=True)
            futures[f] = chunk_paths

        for f in as_completed(futures):
            chunk_paths = futures[f]
            f.result()
            for path in chunk_paths:
                os.unlink(path)
    return
