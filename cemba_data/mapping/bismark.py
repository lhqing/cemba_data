"""
Input: fastq_final_result dataframe
Processes: bismark mapping R1 and R2.
Output: bismark_result dataframe
"""

import pathlib
import pandas as pd
import subprocess
import multiprocessing
import shlex
import functools
import operator


def _parse_bismark_report(report_path_list):
    term_dict = {'Sequences analysed in total': 'total_reads',
                 'Number of alignments with a unique best hit from the different alignments': 'unique_map',
                 'Mapping efficiency': 'mapping_rate',
                 'Sequences with no alignments under any condition': 'unmap',
                 'Sequences did not map uniquely': 'ununique_map',
                 'Sequences which were discarded because genomic sequence could not be extracted': 'no_genome',
                 'CT/CT': 'OT', 'CT/GA': 'OB', 'GA/CT': 'CTOT', 'GA/GA': 'CTOB',
                 'Number of alignments to (merely theoretical) complementary strands being rejected in total': 'reject',
                 "Total number of C's analysed": 'total_c',
                 "Total methylated C's in CpG context": 'total_mcg',
                 "Total methylated C's in CHG context": 'total_mchg',
                 "Total methylated C's in CHH context": 'total_mchh',
                 "Total methylated C's in Unknown context": 'total_mcn',
                 "Total unmethylated C's in CpG context": 'total_cg',
                 "Total unmethylated C's in CHG context": 'total_chg',
                 "Total unmethylated C's in CHH context": 'total_chh',
                 "Total unmethylated C's in Unknown context": 'total_cn',
                 'C methylated in CpG context': 'total_mcg_rate',
                 'C methylated in CHG context': 'total_mchg_rate',
                 'C methylated in CHH context': 'total_mchh_rate',
                 'C methylated in Unknown context (CN or CHN)': 'total_mcn_rate'}

    total_list = []
    for report in report_path_list:
        uid, index_name, read_type = report.name.split('.')[0].split('_')
        with report.open() as rep:
            report_dict = {}
            for line in rep:
                try:
                    start, rest = line.split(':')
                except ValueError:
                    continue  # more or less than 2 after split
                try:
                    report_dict[term_dict[start]] = rest.strip().split('\t')[0].strip('%')
                except KeyError:
                    pass
        report_dict['uid'] = uid
        report_dict['index_name'] = index_name
        report_dict['read_type'] = read_type
        total_list.append(pd.Series(report_dict))
    total_result = pd.DataFrame(total_list)
    return total_result


def bismark(fastq_final_result, out_dir, config):
    bismark_reference = config['bismark']['bismark_reference']
    cores = int(config['bismark']['cores'])
    read_min = int(config['bismark']['read_min'])
    read_max = int(config['bismark']['read_max'])

    bismark_run = functools.partial(subprocess.run,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    encoding='utf8',
                                    check=True)
    # sort by total reads, map large sample first
    sample_dict = {}
    for (uid, index_name), sub_df in fastq_final_result.sort_values('out_reads').groupby(['uid', 'index_name']):
        sample_dict[(uid, index_name)] = sub_df['out_reads'].astype(int).sum()
    sorted_sample = sorted(sample_dict.items(), key=operator.itemgetter(1), reverse=True)

    ran_samples = []
    r1_results = []
    r2_results = []
    pool = multiprocessing.Pool(cores)
    for (uid, index_name), total_reads in sorted_sample:
        if total_reads < read_min:
            print("Drop cell due to too less reads:", uid, index_name, total_reads)
            continue
        if total_reads > read_max:
            print("Drop cell due to too many reads:", uid, index_name, total_reads)
            continue
        ran_samples.append((uid, index_name))
        r1_fastq = pathlib.Path(out_dir) / f'{uid}_{index_name}_R1.trimed.fq.gz'
        r1_cmd = f'bismark {bismark_reference} --bowtie2 {r1_fastq} --pbat -o {out_dir} --temp_dir {out_dir}'
        # each bismark job actually use 250%
        result = pool.apply_async(bismark_run, (shlex.split(r1_cmd),))
        r1_results.append(result)

        r2_fastq = pathlib.Path(out_dir) / f'{uid}_{index_name}_R2.trimed.fq.gz'
        r2_cmd = f'bismark {bismark_reference} --bowtie2 {r2_fastq} -o {out_dir} --temp_dir {out_dir}'
        result = pool.apply_async(bismark_run, (shlex.split(r2_cmd),))
        r2_results.append(result)
    pool.close()
    pool.join()

    if len(ran_samples) != 0:
        report_path_list = []
        for (uid, index_name) in ran_samples:
            report_path_list += list(pathlib.Path(out_dir).glob(f'{uid}_{index_name}*bismark_bt2*report.txt'))
        bismark_result_df = _parse_bismark_report(report_path_list)
        return bismark_result_df
    else:
        return pd.DataFrame()

