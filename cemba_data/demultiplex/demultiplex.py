"""
Demultiplex pipeline
"""

import locale
import logging
import pathlib
import re
import subprocess

import pandas as pd

import cemba_data
from .fastq_dataframe import make_fastq_dataframe
from ..mapping.pipelines import make_snakefile, prepare_run, validate_mapping_config
from ..utilities import snakemake, get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def _demultiplex(fastq_pattern, output_dir, barcode_version, cpu):
    """
    Input raw FASTQ file pattern
    1. automatically parse the name to generate fastq dataframe
    2. use the FASTQ dataframe to generate dir structure for each uid
    3. generate snakefile for demultiplexing each uid
    4. generate final snakefile in output_dir
    5. execute snakefile

    Parameters
    ----------
    fastq_pattern
    output_dir
    barcode_version
    cpu

    Returns
    -------

    """
    output_dir = pathlib.Path(output_dir).absolute()

    # make fastq dataframe
    fastq_df = make_fastq_dataframe(fastq_pattern,
                                    barcode_version=barcode_version,
                                    output_path=output_dir / 'stats' / 'fastq_dataframe.csv')

    # prepare UID sub dir
    snakefile_list = []
    total_stats_list = []
    rule_count = 0
    for uid, uid_df in fastq_df.groupby('uid'):
        # determine index file path
        if barcode_version == 'V1':
            random_index_fasta_path = str(PACKAGE_DIR /
                                          'files/random_index_v1.fa')
        elif barcode_version == 'V2':
            multiplex_group = uid.split('-')[-2]
            random_index_fasta_path = str(
                PACKAGE_DIR / 'files/random_index_v2/'
                              f'random_index_v2.multiplex_group_{multiplex_group}.fa')
        else:
            raise ValueError(f'Got unknown barcode version {barcode_version}.')

        # create a directory for each uid, within this UID, do multiplex and lane merge
        uid_output_dir = output_dir / uid
        uid_output_dir.mkdir(exist_ok=True)
        lane_files_dir = uid_output_dir / 'lanes'
        lane_files_dir.mkdir(exist_ok=True)

        # standardize input fastq name for easier parsing
        raw_dir = uid_output_dir / 'raw'
        raw_dir.mkdir(exist_ok=True)
        for _, row in uid_df.iterrows():
            uid, read_type, lane, old_path = row[[
                'uid', 'read_type', 'lane', 'fastq_path'
            ]]
            new_path = raw_dir / f'{uid}+{lane}+{read_type}.fq.gz'
            subprocess.run(['ln', '-s', old_path, new_path], check=True)
        lanes = list(uid_df['lane'].unique())
        name_str = '{{name}}'

        # make snakefile
        stats_out_list = [
            f'{lane_files_dir}/{uid}-{lane}.demultiplex.stats.txt'
            for lane in lanes
        ]
        total_stats_list += stats_out_list
        rules = ""
        for lane in lanes:
            snake_file_template = f"""
rule demultiplex_{rule_count}:
    input:
        r1_in = f'{raw_dir}/{uid}+{lane}+R1.fq.gz',
        r2_in = f'{raw_dir}/{uid}+{lane}+R2.fq.gz'
    params:
        # Note that you have to use a function to deactivate automatic wildcard expansion 
        # in params strings, e.g., `lambda wildcards: ...`.
        # here the r1/2_out have to have the name_str
        r1_out = lambda wildcards: f'{lane_files_dir}/{uid}-{lane}-{name_str}-R1.fq.gz',
        r2_out = lambda wildcards: f'{lane_files_dir}/{uid}-{lane}-{name_str}-R2.fq.gz'
    output:
        stats_out = '{lane_files_dir}/{uid}-{lane}.demultiplex.stats.txt'
    shell:
        "cutadapt -Z -e 0.01 --no-indels -g file:{random_index_fasta_path} "
        "-o {{params.r1_out}} -p {{params.r2_out}} {{input.r1_in}} {{input.r2_in}} > {{output.stats_out}}"
    """
            rule_count += 1
            rules += snake_file_template

        snake_file_path = lane_files_dir / 'Snakefile'
        with open(snake_file_path, 'w') as f:
            f.write(rules)
        snakefile_list.append(f'{uid}/lanes/Snakefile')

    # make final snakefile for demultiplex step
    final_rules = ''
    for path in snakefile_list:
        final_rules += f'include: "{path}"\n'
    # final rules
    final_rules += f"""
rule final:
    input: {total_stats_list}
"""
    final_snake_path = output_dir / 'Snakefile_demultiplex'
    with open(final_snake_path, 'w') as f:
        f.write(final_rules)

    print('Demultiplexing raw FASTQ')
    snakemake(workdir=output_dir, snakefile=final_snake_path, cores=cpu)
    return


def _merge_lane(output_dir, cpu):
    output_dir = pathlib.Path(output_dir).absolute()
    fastq_df = pd.read_csv(output_dir / 'stats' / 'fastq_dataframe.csv')
    snakefile_list = []
    total_output_list = []
    rule_uid = 0
    # prepare snakefile in each uid
    for uid in fastq_df['uid'].unique():
        uid_output_dir = output_dir / uid
        lanes_dir = uid_output_dir / 'lanes'
        fastq_dir = uid_output_dir / 'fastq'
        fastq_dir.mkdir(exist_ok=True)

        # prepare demultiplex results cell_fastq_df
        records = []
        for path in lanes_dir.glob('*fq.gz'):
            *uid, lane, index_name, read_type = path.name[:-6].split('-')
            uid = '-'.join(uid)
            cell_id = f'{uid}-{index_name}'
            records.append([cell_id, lane, read_type, str(path)])
        cell_fastq_df = pd.DataFrame(
            records, columns=['cell_id', 'index_name', 'read_type', 'fastq_path'])

        # prepare snakefile for each cell_id * read_type
        rules = ''
        output_paths = []
        for (cell_id, read_type), sub_df in cell_fastq_df.groupby(['cell_id', 'read_type']):
            input_paths = list(sub_df['fastq_path'])
            output_path = fastq_dir / f'{cell_id}-{read_type}.fq.gz'

            snake_file_template = f"""
rule merge_{rule_uid}:
    input: 
        {input_paths}
    output: 
        "{output_path}"
    shell:
        "gzip -cd {{input}} | gzip -5 > {{output}} && rm -f {{input}}"

"""
            rule_uid += 1
            rules += snake_file_template
            output_paths.append(str(output_path))

        snakefile_path = uid_output_dir / 'Snakefile'
        with open(snakefile_path, 'w') as f:
            f.write(rules)
        snakefile_list.append(snakefile_path)
        total_output_list += output_paths

    # prepare final snakefile
    final_rules = ''
    for path in snakefile_list:
        final_rules += f'include: "{path}"\n'
    # final rules
    final_rules += f"""
rule final:
    input: {total_output_list}
"""
    final_snake_path = output_dir / 'Snakefile_merge_lane'
    with open(final_snake_path, 'w') as f:
        f.write(final_rules)

    print('Merging lanes to get cell FASTQ')
    snakemake(workdir=output_dir, snakefile=final_snake_path, cores=cpu)
    return


def _parse_index_fasta(fasta_path):
    records = {}
    with open(fasta_path) as f:
        key_line = True
        for line in f:
            if key_line:
                key = line.lstrip('>').rstrip('\n')
                key_line = False
            else:
                value = line.lstrip('^').rstrip('\n')
                records[key] = value
                key_line = True
    return records


def _read_cutadapt_result(stat_path):
    """
    Parser of cutadapt output
    """
    with open(stat_path) as f:
        p = re.compile(
            r"Sequence: .+; Type: .+; Length: \d+; Trimmed: \d+ times")
        series = []
        total_pairs = -1
        for line in f:
            if line.startswith('Total read pairs processed'):
                # some weird transform of cutadapt outputs...
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
        total_df = pd.DataFrame(series)
        total_df['Trimmed'] = total_df['Trimmed'].apply(
            lambda c: c.split(' ')[0]).astype(int)
        total_df['TotalPair'] = total_pairs
        total_df['Ratio'] = total_df['Trimmed'] / total_pairs
    return total_df


def _summarize_demultiplex(output_dir, barcode_version):
    output_dir = pathlib.Path(output_dir).absolute()
    output_path = output_dir / 'stats' / 'demultiplex.stats.csv'
    barcode_version = barcode_version.upper()

    # get index info
    if barcode_version == 'V1':
        random_index_fasta_path = str(PACKAGE_DIR /
                                      'files/random_index_v1.fa')
    elif barcode_version == 'V2':
        # here we don't need to worry about the multiplex_group issue,
        # because we just need a index_name to index_seq map
        # we've considered this during demultiplex
        random_index_fasta_path = str(
            PACKAGE_DIR / 'files/random_index_v2/random_index_v2.fa')
    else:
        raise ValueError(
            f'Unknown version name {barcode_version} in multiplexIndex section of the config file.'
        )
    index_seq_dict = _parse_index_fasta(random_index_fasta_path)
    index_name_dict = {v: k for k, v in index_seq_dict.items()}

    # read the demultiplex stats, its per lane, so need to sum up lane together of each uid and index name
    # but R1 R2 is demultiplexed together, so this table don't separate R1 R2
    stat_list = []
    stat_path_list = list(output_dir.glob('*/lanes/*demultiplex.stats.txt'))
    for path in stat_path_list:
        single_df = _read_cutadapt_result(path)
        *uid, suffix = path.name.split('-')
        lane = suffix.split('.')[0]
        uid = '-'.join(uid)
        single_df['uid'] = uid
        single_df['lane'] = lane
        single_df['index_name'] = single_df['Sequence'].map(index_name_dict)
        assert single_df['index_name'].isna().sum() == 0
        stat_list.append(single_df)
    total_demultiplex_stats = pd.concat(stat_list)

    # calculate cell level table
    total_demultiplex_stats['cell_id'] = total_demultiplex_stats[
                                             'uid'] + '-' + total_demultiplex_stats['index_name']

    cell_table = total_demultiplex_stats.groupby('cell_id').agg({
        'Trimmed': 'sum',
        'TotalPair': 'sum',
        'index_name': lambda i: i.unique()[0],
        'uid': lambda i: i.unique()[0]
    })
    cell_table.rename(columns={
        'Trimmed': 'CellInputReadPairs',
        'TotalPair': 'MultiplexedTotalReadPairs',
        'index_name': 'IndexName',
        'uid': 'UID'
    },
        inplace=True)
    cell_table['CellBarcodeRate'] = cell_table[
                                        'CellInputReadPairs'] / cell_table['MultiplexedTotalReadPairs']
    cell_table['BarcodeVersion'] = barcode_version
    cell_table.to_csv(output_path)
    return


def _final_cleaning(output_dir):
    """
    remove intermediate files from demultiplex
    """
    output_dir = pathlib.Path(output_dir)

    delete_patterns = [f'Snakefile_*', '*/lanes', '*/raw', '*/Snakefile', '*/fastq/*-unknown-R*.fq.gz', '.snakemake']

    total_paths = []
    for pattern in delete_patterns:
        total_paths += list(map(str, output_dir.glob(pattern)))

    subprocess.run(['rm', '-rf'] + total_paths, check=True)
    return


def _skip_abnormal_fastq_pairs(output_dir):
    demultiplex_df = pd.read_csv(output_dir / 'stats/demultiplex.stats.csv', index_col=0)
    config = get_configuration(output_dir / 'mapping_config.ini')
    total_read_pairs_min = int(config['total_read_pairs_min'])
    total_read_pairs_max = int(config['total_read_pairs_max'])

    too_large = demultiplex_df['CellInputReadPairs'] > total_read_pairs_max
    too_small = demultiplex_df['CellInputReadPairs'] < total_read_pairs_min
    judge = too_small | too_large
    unmapped_cells = demultiplex_df[judge]
    print(f'Skip {too_small.sum()} cells due to too less input read pairs (< {total_read_pairs_min})')
    print(f'Skip {too_large.sum()} cells due to too large input read pairs (> {total_read_pairs_max})')

    for cell_id, row in unmapped_cells.iterrows():
        uid = row['UID']
        skipped_dir = output_dir / uid / 'fastq/skipped/'
        skipped_dir.mkdir(exist_ok=True)

        # move both R1 R2 to skipped files, it will not be included in Snakefile
        for read_type in ['R1', 'R2']:
            fastq_path = output_dir / uid / f'fastq/{cell_id}-{read_type}.fq.gz'
            new_path = skipped_dir / f'{cell_id}-{read_type}.fq.gz'
            # if CellInputReadPairs = 0, the FASTQ file do not actually exist, but it does have a row in metadata.
            if fastq_path.exists():
                subprocess.run(['mv', str(fastq_path), str(new_path)], check=True)

    # save UID total input reads, for command order
    uid_order = demultiplex_df[~judge].groupby(
        'UID')['CellInputReadPairs'].sum().sort_values(
        ascending=False)
    uid_order.to_csv(output_dir / 'stats/UIDTotalCellInputReadPairs.csv', header=False)
    return


SUPPORTED_TECHNOLOGY = ['mc', 'mct', 'm3c']


def demultiplex_pipeline(fastq_pattern, output_dir, config_path, cpu):
    cpu = int(cpu)
    merge_cpu = min(48, cpu)
    demultiplex_cpu = min(32, cpu)

    output_dir = pathlib.Path(output_dir).absolute()
    if output_dir.exists():
        raise FileExistsError('output_dir already exists, to prevent conflicts, '
                              'use another output_dir or delete the existing output_dir first.')
    else:
        output_dir.mkdir(parents=True)
        (output_dir / 'stats').mkdir()

    config = get_configuration(config_path)
    new_config_path = output_dir / 'mapping_config.ini'
    subprocess.run(f'cp {config_path} {new_config_path}', shell=True, check=True)
    barcode_version = config['barcode_version']
    # validate config file first before demultiplex
    validate_mapping_config(output_dir)

    _demultiplex(
        fastq_pattern=fastq_pattern,
        output_dir=output_dir,
        barcode_version=barcode_version,
        cpu=demultiplex_cpu)
    _merge_lane(output_dir=output_dir, cpu=merge_cpu)
    _summarize_demultiplex(output_dir=output_dir, barcode_version=barcode_version)
    _final_cleaning(output_dir=output_dir)
    _skip_abnormal_fastq_pairs(output_dir=output_dir)
    make_snakefile(output_dir=output_dir)

    # this is just a convenient step, so I fix the parameters here
    # users should change the resulting batch submission
    # or generate by themselves if they want different setting.
    prepare_run(output_dir)
    return
