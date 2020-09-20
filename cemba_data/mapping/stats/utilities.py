import pathlib

import pandas as pd

from ...utilities import parse_mc_pattern


def parse_trim_fastq_stats(stat_path):
    # example trim fastq stats
    """
status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
0	OK	1490	213724	0	0	0	1490	4	0	213712
1	status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
2	OK	1490	213712	0	0	0	1482	0	1300	182546
"""
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    trim_stats = pd.read_csv(stat_path, sep='\t')
    trim_stats = trim_stats.iloc[[0, 2], :].reset_index()  # skip the duplicated title row

    data = pd.Series({
        f'{read_type}InputReads': trim_stats['in_reads'][0],
        f'{read_type}InputReadsBP': trim_stats['in_bp'][0],
        f'{read_type}WithAdapters': trim_stats['w/adapters'][0],
        f'{read_type}QualTrimBP': trim_stats['qualtrim_bp'][1],
        f'{read_type}TrimmedReads': trim_stats['out_reads'][1],
        f'{read_type}TrimmedReadsBP': trim_stats['out_bp'][1],
        f'{read_type}TrimmedReadsRate': int(trim_stats['out_reads'][1]) / int(trim_stats['in_reads'][0])
    }, name=cell_id)
    return data


def parse_bismark_report(stat_path):
    """
    parse bismark output
    """
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    term_dict = {
        'Number of alignments with a unique best hit from the different alignments': f'{read_type}UniqueMappedReads',
        'Mapping efficiency': f'{read_type}MappingRate',
        'Sequences with no alignments under any condition': f'{read_type}UnmappedReads',
        'Sequences did not map uniquely': f'{read_type}UnuniqueMappedReads',
        'CT/CT': f'{read_type}OT',
        'CT/GA': f'{read_type}OB',
        'GA/CT': f'{read_type}CTOT',
        'GA/GA': f'{read_type}CTOB',
        "Total number of C's analysed": f'{read_type}TotalC',
        'C methylated in CpG context': f'{read_type}TotalmCGRate',
        'C methylated in CHG context': f'{read_type}TotalmCHGRate',
        'C methylated in CHH context': f'{read_type}TotalmCHHRate'}

    with open(stat_path) as rep:
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
    return pd.Series(report_dict, name=cell_id)


def parse_deduplicate_stat(stat_path):
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    try:
        dedup_result_series = pd.read_csv(stat_path, comment='#', sep='\t').T[0]
        rename_dict = {
            'UNPAIRED_READS_EXAMINED': f'{read_type}MAPQFilteredReads',
            'UNPAIRED_READ_DUPLICATES': f'{read_type}DuplicatedReads',
            'PERCENT_DUPLICATION': f'{read_type}DuplicationRate'
        }
        dedup_result_series = dedup_result_series.loc[rename_dict.keys()].rename(rename_dict)

        dedup_result_series[f'{read_type}FinalBismarkReads'] = dedup_result_series[f'{read_type}MAPQFilteredReads'] - \
                                                               dedup_result_series[f'{read_type}DuplicatedReads']
        dedup_result_series.name = cell_id
    except pd.errors.EmptyDataError:
        # if a BAM file is empty, picard matrix is also empty
        dedup_result_series = pd.Series({f'{read_type}MAPQFilteredReads': 0,
                                         f'{read_type}DuplicatedReads': 0,
                                         f'{read_type}FinalBismarkReads': 0,
                                         f'{read_type}DuplicationRate': 0}, name=cell_id)
    return dedup_result_series


def generate_allc_stats(output_dir, config):
    output_dir = pathlib.Path(output_dir).absolute()
    allc_stats_dict = {p.name.split('.')[0]: p for p in output_dir.glob('allc/*count.csv')}

    patterns = config['mc_stat_feature'].split(' ')
    patterns_alias = config['mc_stat_alias'].split(' ')
    pattern_translate = {k: v for k, v in zip(patterns, patterns_alias)}

    # real all cell stats
    total_stats = []
    for cell_id, path in allc_stats_dict.items():
        allc_stat = pd.read_csv(path, index_col=0)
        allc_stat['cell_id'] = cell_id
        total_stats.append(allc_stat)
    total_stats = pd.concat(total_stats)
    cell_genome_cov = pd.Series(total_stats.set_index('cell_id')['genome_cov'].to_dict())
    # aggregate into patterns
    cell_records = []
    for pattern in pattern_translate.keys():
        contexts = parse_mc_pattern(pattern)
        pattern_stats = total_stats[total_stats.index.isin(contexts)]
        cell_level_data = pattern_stats.groupby('cell_id')[['mc', 'cov']].sum()
        cell_level_data['frac'] = cell_level_data['mc'] / cell_level_data['cov']

        # prettify col name
        _pattern = pattern_translate[pattern]
        cell_level_data = cell_level_data.rename(
            columns={'frac': f'{_pattern}Frac',
                     'mc': f'{_pattern}mC',
                     'cov': f'{_pattern}Cov'})
        cell_records.append(cell_level_data)
    final_df = pd.concat(cell_records, axis=1, sort=True).reindex(allc_stats_dict.keys())
    final_df['GenomeCov'] = cell_genome_cov
    final_df.index.name = 'cell_id'
    return final_df
