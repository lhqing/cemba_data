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


def star_mapping(input_dir, output_dir, config):
    """
    STAR RNA reads mapping
    """
    output_dir = pathlib.Path(output_dir)
    input_dir = pathlib.Path(input_dir)
    fastq_qc_records = pd.read_csv(input_dir / 'fastq_qc.records.csv',
                                   index_col=['uid', 'index_name', 'read_type'], squeeze=True)
    if isinstance(config, str):
        config = get_configuration(config)

    star_reference = config['star']['star_reference']
    read_min = int(config['star']['read_min'])
    read_max = int(config['star']['read_max'])

    # Do not allow to change threads
    # threads = int(config['star']['threads'])
    threads = 6

    try:
        fastq_qc_stats_path = input_dir / 'fastq_qc.stats.csv'
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

        # for RNA part, only map R1
        # TODO Do pair end mapping...
        r1_fastq = fastq_qc_records[(uid, index_name, 'R1')]
        output_prefix = output_dir / (pathlib.Path(r1_fastq).name[:-6] + '.STAR.')
        output_path = output_dir / (pathlib.Path(r1_fastq).name[:-6] + '.STAR.Aligned.out.bam')  # STAR convention
        star_cmd = f'STAR --runThreadN {threads} ' \
                   f'--genomeDir {star_reference} ' \
                   f'--alignEndsType EndToEnd ' \
                   f'--genomeLoad LoadAndKeep ' \
                   f'--outSAMstrandField intronMotif ' \
                   f'--outSAMtype BAM Unsorted ' \
                   f'--outSAMunmapped Within ' \
                   f'--outSAMattributes NH HI AS NM MD ' \
                   f'--sjdbOverhang 100 ' \
                   f'--outFilterType BySJout ' \
                   f'--outFilterMultimapNmax 20 ' \
                   f'--alignSJoverhangMin 8 ' \
                   f'--alignSJDBoverhangMin 1 ' \
                   f'--outFilterMismatchNmax 999 ' \
                   f'--outFilterMismatchNoverLmax 0.04 ' \
                   f'--alignIntronMin 20 ' \
                   f'--alignIntronMax 1000000 ' \
                   f'--alignMatesGapMax 1000000 ' \
                   f'--outFileNamePrefix {output_prefix} ' \
                   f'--readFilesIn {r1_fastq} ' \
                   f'--readFilesCommand gzip -cd'
        records.append([uid, index_name, output_path])
        command_list.append(star_cmd)

    fold_command_list = []
    # fold STAR cmd and add genome load genome remove
    command_per_script = 50
    genome_load_cmd = f'STAR --genomeDir {star_reference} --genomeLoad LoadAndExit'
    genome_remove_cmd = f'STAR --genomeDir {star_reference} --genomeLoad Remove'
    for i in range(0, len(command_list), command_per_script):
        commands_str = '\n'.join(command_list[i: i+command_per_script])
        total_command = f'{genome_load_cmd}\n{commands_str}\n{genome_remove_cmd}'
        fold_command_list.append(total_command)

    with open(output_dir / 'star_mapping.command.txt', 'w') as f:
        f.write('\n'.join(fold_command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'bam_path'])
    record_df.to_csv(output_dir / 'star_mapping.records.csv', index=None)
    return record_df, fold_command_list


def _parse_star_log(log_path):
    with open(log_path) as f:
        record = {}
        for line in f:
            if '|' not in line:
                continue
            name, value = line.split('|')
            name = name.strip(' ')
            value = value.strip(' \t\n')
            record[name] = value
    return pd.Series(record)


def summarize_star_mapping(bam_dir):
    bam_dir = pathlib.Path(bam_dir)
    output_path = bam_dir / 'star_mapping.stats.csv'
    if output_path.exists():
        return str(output_path)

    records = []
    star_stat_list = list(bam_dir.glob('*.Log.final.out'))
    for path in star_stat_list:
        report_series = _parse_star_log(path)

        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        report_series['uid'] = uid
        report_series['index_name'] = index_name
        report_series['read_type'] = read_type
        records.append(report_series)
        path = str(path)
        subprocess.run(['rm', '-f',
                        path,
                        path[:-10] + '.out',
                        path[:-10] + '.progress.out',
                        path[:-14] + '.SJ.out.tab'])
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path)
    return str(output_path)
