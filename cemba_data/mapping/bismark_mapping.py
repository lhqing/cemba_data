import logging
import operator
import pathlib
import subprocess

import pandas as pd

import cemba_data
from .utilities import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def bismark_mapping(input_dir, output_dir, config):
    """
    reads level QC and trimming. R1 R2 separately.
    """
    output_dir = pathlib.Path(output_dir)
    input_dir = pathlib.Path(input_dir)
    fastq_qc_records = pd.read_csv(input_dir / 'fastq_qc.records.csv',
                                   index_col=['uid', 'index_name', 'read_type'], squeeze=True)
    if isinstance(config, str):
        config = get_configuration(config)

    bismark_reference = config['bismark']['bismark_reference']
    read_min = int(config['bismark']['read_min'])
    read_max = int(config['bismark']['read_max'])
    try:
        remove_fastq_input = config['bismark']['remove_fastq_input']
        if 'mct' in config and (config['mct'].lower() in ['true', 't', 'yes', 'y']):
            # for mCT seq, do not remove filtered fastq, cause STAR mapping use it
            remove_fastq_input = False

    except KeyError:
        remove_fastq_input = 'True'
        log.warning('remove_fastq_input not found in config.ini file, you are using the old version'
                    'please update your config.ini use yap default-mapping-config')

    fastq_qc_stats_path = input_dir / 'fastq_qc.stats.csv'
    try:
        fastq_qc_stats = pd.read_csv(fastq_qc_stats_path, index_col=0)
        # sort by total reads, map large sample first
        sample_dict = {}
        for (uid, index_name), sub_df in fastq_qc_stats.sort_values('out_reads').groupby(['uid', 'index_name']):
            sample_dict[(uid, index_name)] = sub_df['out_reads'].astype(int).sum()
        sorted_sample = sorted(sample_dict.items(), key=operator.itemgetter(1), reverse=True)
    except FileNotFoundError:
        # not stats means in dry run mode, will generate command for every record
        sorted_sample = []
        for (uid, index_name), _ in fastq_qc_records.groupby(['uid', 'index_name']):
            sorted_sample.append([(uid, index_name), -1])

    records = []
    command_list = []
    for (uid, index_name), total_reads in sorted_sample:
        if index_name == 'unknown':
            continue
        elif total_reads == -1:
            pass
        elif total_reads < read_min:
            log.info(f"Drop cell due to too less reads: {uid} {index_name}, total reads {total_reads}")
            continue
        elif total_reads > read_max:
            log.info(f"Drop cell due to too many reads: {uid} {index_name}, total reads {total_reads}")
            continue
        r1_fastq = fastq_qc_records[(uid, index_name, 'R1')]
        remove_cmd = f' && rm -f {r1_fastq}' if (remove_fastq_input.lower() in ['y', 'yes', 't', 'true']) else ''
        # run R1 with pbat mode
        r1_cmd = f'bismark {bismark_reference} --bowtie2 {r1_fastq} ' \
                 f'--pbat -o {output_dir} --temp_dir {output_dir}{remove_cmd}'
        r1_output_path = output_dir / (pathlib.Path(r1_fastq).name[:-6] + '_bismark_bt2.bam')  # bismark convention
        command_list.append(r1_cmd)
        records.append([uid, index_name, 'R1', r1_output_path])

        r2_fastq = fastq_qc_records[(uid, index_name, 'R2')]
        # run R2 with normal mode
        remove_cmd = f' && rm -f {r2_fastq}' if (remove_fastq_input.lower() in ['y', 'yes', 't', 'true']) else ''
        r2_cmd = f'bismark {bismark_reference} --bowtie2 {r2_fastq} ' \
                 f'-o {output_dir} --temp_dir {output_dir}{remove_cmd}'
        r2_output_path = output_dir / (pathlib.Path(r2_fastq).name[:-6] + '_bismark_bt2.bam')  # bismark convention
        command_list.append(r2_cmd)
        records.append([uid, index_name, 'R2', r2_output_path])

    with open(output_dir / 'bismark_mapping.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'bam_path'])
    record_df.to_csv(output_dir / 'bismark_mapping.records.csv', index=None)
    return record_df, command_list


def _parse_bismark_report(report_path):
    """
    Some ugly parser for bismark_mapping report... Hope Bismark won't change...
    # TODO make this independent to bismark_mapping
    """
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

    with report_path.open() as rep:
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
    return pd.Series(report_dict)


def summarize_bismark_mapping(output_dir):
    bam_dir = pathlib.Path(output_dir)
    output_path = bam_dir / 'bismark_mapping.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    bismark_stat_list = list(bam_dir.glob('*_bt2_SE_report.txt'))
    for path in bismark_stat_list:
        report_series = _parse_bismark_report(path)

        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        report_series['uid'] = uid
        report_series['index_name'] = index_name
        report_series['read_type'] = read_type
        records.append(report_series)
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path, index=None)
    return str(output_path)
