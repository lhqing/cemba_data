"""
Input: bismark_result dataframe
Processes: sort, dedup and mapq filter on bam; merge R1 R2 bam.
Output: bam_result dataframe
"""

import pathlib
import pandas as pd
import subprocess
import multiprocessing
import shlex
import glob


def _process_bam(cmd_list):
    return [subprocess.run(shlex.split(cmd),
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           encoding='utf8')
            for cmd in cmd_list]


def bam_qc(bismark_result, out_dir, config):
    cores = int(config['bamFilter']['cores'])
    mapq_threshold = config['bamFilter']['mapq_threshold']

    # process bam
    pool = multiprocessing.Pool(cores)
    results = []
    for i, line in bismark_result.iterrows():
        uid, index_name, read_type = line[['uid', 'index_name', 'read_type']]
        # file path
        bismark_bam = str(
            list(pathlib.Path(out_dir).glob(f'{uid}_{index_name}_{read_type}*_bismark_bt2.bam'))[0].absolute())
        sort_bam = bismark_bam[:-3] + 'sort.bam'
        dedup_bam = bismark_bam[:-3] + 'dedup.bam'
        dedup_matrix = bismark_bam[:-3] + 'dedup.matrix.txt'
        filter_bam = bismark_bam[:-3] + 'filter.bam'
        # command
        sort_cmd = f'samtools sort -o {sort_bam} --threads 2 {bismark_bam}'
        dedup_cmd = f'picard MarkDuplicates I={sort_bam} O={dedup_bam} M={dedup_matrix} REMOVE_DUPLICATES=true'
        filter_cmd = f'samtools view -b -h -q {mapq_threshold} -o {filter_bam} {dedup_bam}'
        cmd_list = [sort_cmd, dedup_cmd, filter_cmd]
        result = pool.apply_async(_process_bam, (cmd_list,))
        results.append(result)
    pool.close()
    pool.join()

    # get qc stats
    qc_result = []
    for i, line in bismark_result.iterrows():
        uid, index_name, read_type = line[['uid', 'index_name', 'read_type']]
        dedup_matrix = list(pathlib.Path(out_dir).glob(f'{uid}_{index_name}_{read_type}*dedup.matrix.txt'))[0]
        filter_bam = list(pathlib.Path(out_dir).glob(f'{uid}_{index_name}_{read_type}*filter.bam'))[0]
        count_filter_bam_cmd = f'samtools view -c {filter_bam}'
        count_result = subprocess.run(shlex.split(count_filter_bam_cmd),
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        remain_reads = int(count_result.stdout.strip())
        header = True
        lines = []
        for _line in dedup_matrix.open().read().split('\n'):
            if _line.startswith('LIBRARY'):
                header = False
            if header:
                continue
            ll = _line.strip().split('\t')
            if len(ll) != 1:
                lines.append(ll)
        s = pd.Series({k: v for k, v in zip(*lines)})
        s['uid'] = uid
        s['index_name'] = index_name
        s['read_type'] = read_type
        s['out_reads'] = remain_reads
        qc_result.append(s)
    bam_result_df = pd.DataFrame(qc_result)

    for (uid, index_name), sub_df in bam_result_df.groupby(['uid', 'index_name']):
        # merge R1 R2 bam
        merge_bam = pathlib.Path(out_dir) / f'{uid}_{index_name}.final.bam'
        bam_list = ' '.join(glob.glob(f'{out_dir}/{uid}_{index_name}_R*filter.bam'))
        merge_cmd = f'samtools merge -f {merge_bam} {bam_list}'
        subprocess.run(shlex.split(merge_cmd),
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
        # clean up
        # only keep filtered bam and bismark report
        remove_file_list = [str(p) for p in pathlib.Path(out_dir).glob(f'{uid}_{index_name}_R*.*bismark*')]
        remove_file_str = ' '.join(remove_file_list)
        rm_cmd = f'rm -f {remove_file_str}'
        subprocess.run(shlex.split(rm_cmd), encoding='utf8')
    return bam_result_df
