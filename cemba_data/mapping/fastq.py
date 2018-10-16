import pathlib
import pandas as pd
import subprocess
import re
import locale
import functools
import multiprocessing
import shlex


def _make_command_dataframe(fastq_dataframe, out_dir, config):
    multiplex_index_dict = config['multiplexIndex']
    overlap = int(config['demultiplex']['overlap'])
    anchor = bool(config['demultiplex']['anchor'])
    adapter_pos = int(config['demultiplex']['adapter_pos'])

    # make adapter parms
    adapter_type = '-g' if int(adapter_pos) == 5 else '-a'
    if anchor:
        multiplex_index_dict = {k: 'X' + v for k, v in multiplex_index_dict.items()}
    adapter_parms = ' '.join([f'{adapter_type} {k}={v}' for k, v in multiplex_index_dict.items()])

    # make cmd_list
    required_cols = ('uid', 'lane', 'read-type', 'fastq-path')
    for col in required_cols:
        if col not in fastq_dataframe.columns:
            raise ValueError(col, 'not in fastq dataframe')
    # standardize read_type
    fastq_dataframe['read-type'] = fastq_dataframe['read-type'].apply(lambda i: 'R2' if '2' in str(i) else 'R1')

    records = []
    for (uid, lane), sub_df in fastq_dataframe.groupby(['uid', 'lane']):
        tmp_sub_df = sub_df.set_index('read-type')
        r1_in = tmp_sub_df.loc['R1', 'fastq-path']
        r2_in = tmp_sub_df.loc['R2', 'fastq-path']
        r1_out = pathlib.Path(out_dir) / (f"{uid}_{lane}" + "_{name}_R1.fq.gz")
        r2_out = pathlib.Path(out_dir) / (f"{uid}_{lane}" + "_{name}_R2.fq.gz")
        cmd = f"cutadapt {adapter_parms} -O {overlap} -o {r1_out.absolute()} -p {r2_out.absolute()} {r1_in} {r2_in}"
        records.append([uid, lane, cmd])

    cmd_df = pd.DataFrame(records, columns=['uid', 'lane', 'cmd'])
    return cmd_df


def _read_cutadapt_result(result):
    p = re.compile(r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times.")
    p_name = re.compile(r'First read: Adapter AD\d+')
    series = []
    adapter_names = []
    for line in result.stdout.split('\n'):
        if line.startswith('Total read pairs processed'):
            locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
            total_pairs = locale.atoi(line.split(' ')[-1])

        m = p.search(line)
        if m is not None:
            result_dict = {}
            for i in m.group().split('; '):
                k, v = i.split(': ')
                result_dict[k] = v
            result_series = pd.Series(result_dict)
            series.append(result_series)
        m = p_name.search(line)
        if m is not None:
            adapter_names.append(m.group().split(' ')[-1])
    total_df = pd.DataFrame(series)
    total_df['Trimmed'] = total_df['Trimmed'].apply(lambda i: i.split(' ')[0]).astype(int)
    total_df['Adepter'] = adapter_names
    total_df.set_index('Adepter', inplace=True)
    total_df['TotalPair'] = total_pairs
    total_df['Ratio'] = total_df['Trimmed'] / total_pairs
    return total_df


def _demultiplex(fastq_dataframe, out_dir, config):
    multiplex_index_dict = config['multiplexIndex']
    multiplex_index_map = {v: k for k, v in multiplex_index_dict.items()}
    cmd_df = _make_command_dataframe(fastq_dataframe, out_dir, config)

    cutadapt_run = functools.partial(subprocess.run,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE,
                                     encoding='utf8')

    total_results = []
    for uid, single_uid_df in cmd_df.groupby('uid'):
        pool = multiprocessing.Pool(cmd_df.shape[0])
        results = []
        for i, row in single_uid_df.iterrows():
            result = pool.apply_async(cutadapt_run, (shlex.split(row['cmd']),))
            results.append(result)
        pool.close()
        pool.join()

        for result, lane in zip(results, single_uid_df['lane']):
            result_df = _read_cutadapt_result(result.get())
            result_df['lane'] = lane
            result_df['uid'] = uid
            result_df['index_name'] = result_df['Sequence'].apply(lambda ind: multiplex_index_map[ind])
            total_results.append(result_df)
    total_result_df = pd.concat(total_results, ignore_index=True)

    # clean small fastq
    # TODO
    return total_result_df


def _fastq_qc(demultiplex_result, out_dir, config):
    pigz_cores = 4
    cutadapt_cores = 4

    r1_adapter = config['fastqTrim']['r1_adapter']
    r2_adapter = config['fastqTrim']['r1_adapter']
    length_threshold = config['fastqTrim']['length_threshold']
    quality_threshold = config['fastqTrim']['quality_threshold']
    left_cut = config['fastqTrim']['left_cut']
    right_cut = config['fastqTrim']['right_cut']
    overlap = config['fastqTrim']['overlap']
    total_reads_threshold = int(config['fastqTrim']['total_reads_threshold'])

    results = []
    for (uid, index_name), sub_df in demultiplex_result.groupby(['uid', 'index_name']):
        sample_demultiplex_total = sub_df['Trimmed'].sum()
        if sample_demultiplex_total < total_reads_threshold:
            continue
        # process R1
        r1_path_pattern = f'{out_dir}/{uid}_L*_{index_name}_R1.fq.gz'
        r1_out = f'{out_dir}/{uid}_{index_name}_R1.trimed.fq.gz'
        r1_cmd = f'pigz -cd -p {pigz_cores} {r1_path_pattern} | ' \
                 f'cutadapt -j {cutadapt_cores} --report=minimal -O {overlap} ' \
                 f'-q {quality_threshold} -u {left_cut} ' \
                 f'-u -{right_cut} -m {length_threshold} ' \
                 f'-a {r1_adapter} -o {r1_out} -'
        r1_result = subprocess.run(r1_cmd, stdout=subprocess.PIPE, encoding='utf8', shell=True)
        # get R1 result stat
        lines = []
        for line in r1_result.stdout.split('\n'):
            l = line.split('\t')
            if len(l) > 1:
                lines.append(l)
        s = pd.Series({name: number for name, number in zip(*lines)})
        s['uid'] = uid
        s['index_name'] = index_name
        s['read_type'] = 'R1'
        results.append(s)

        # process R2
        r2_path_pattern = f'{out_dir}/{uid}_L*_{index_name}_R2.fq.gz'
        r2_out = f'{out_dir}/{uid}_{index_name}_R2.trimed.fq.gz'
        r2_cmd = f'pigz -cd -p {pigz_cores} {r2_path_pattern} | ' \
                 f'cutadapt -j {cutadapt_cores} --report=minimal -O {overlap} ' \
                 f'-q {quality_threshold} -u {left_cut} ' \
                 f'-u -{right_cut} -m {length_threshold} ' \
                 f'-a {r2_adapter} -o {r2_out} -'
        r2_result = subprocess.run(r2_cmd, stdout=subprocess.PIPE, encoding='utf8', shell=True)
        # get R2 result stat
        lines = []
        for line in r2_result.stdout.split('\n'):
            l = line.split('\t')
            if len(l) > 1:
                lines.append(l)
        s = pd.Series({name: number for name, number in zip(*lines)})
        s['uid'] = uid
        s['index_name'] = index_name
        s['read_type'] = 'R2'
        results.append(s)

    fastq_final_result = pd.DataFrame(results)
    fastq_final_result['out_reads_rate'] = \
        fastq_final_result['out_reads'].astype(int) / fastq_final_result['in_reads'].astype(int)
    fastq_final_result['out_bp_rate'] = \
        fastq_final_result['out_reads'].astype(int) / fastq_final_result['in_reads'].astype(int)

    # clean up
    for (uid, index_name), sub_df in demultiplex_result.groupby(['uid', 'index_name']):
        r_path_pattern = f'{out_dir}/{uid}_L*_{index_name}_R*.fq.gz'
        r_rm_cmd = f'rm -f {r_path_pattern}'
        subprocess.run(r_rm_cmd, shell=True)

    return fastq_final_result
