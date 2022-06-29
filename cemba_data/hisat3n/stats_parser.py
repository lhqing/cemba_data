import pathlib

import pandas as pd


def parse_hisat_report(stat_path):
    """
    parse hisat3n output
    """
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    term_dict = {
        'Total pairs': f'TotalPairs',
        'Aligned concordantly 1 time': f'UniqueMappedPairs',
        'Aligned discordantly 1 time': f'discordantlyUniqueMappedPairs',
        'Aligned 1 time': f'UniqueMappedReads',
    }

    with open(stat_path) as rep:
        report_dict = {}
        for line in rep:
            try:
                start, rest = line.split(':')
                start = start.strip()
            except ValueError:
                continue  # more or less than 2 after split
            try:
                report_dict[term_dict[start]] = rest.strip().split(' ')[0].strip('%')
            except KeyError:
                pass
        mappingrate = (int(report_dict[f'UniqueMappedPairs']) * 2 + int(
            report_dict[f'discordantlyUniqueMappedPairs']) * 2 + int(report_dict[f'UniqueMappedReads'])) / (
                                  int(report_dict[f'TotalPairs']) * 2)
        report_dict[f'MappingRate'] = round(mappingrate, 2)
    return pd.Series(report_dict, name=cell_id)


def parse_deduplicate_stat(stat_path):
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    try:
        dedup_result_series = pd.read_csv(stat_path, comment='#', sep='\t').T[0]
        rename_dict = {
            'UNPAIRED_READS_EXAMINED': f'FilteredReads',
            'READ_PAIRS_EXAMINED': f'FilteredPairs',
            'UNPAIRED_READ_DUPLICATES': f'DuplicatedReads',
            'READ_PAIR_DUPLICATES': f'DuplicatedPairs',
            'PERCENT_DUPLICATION': f'DuplicationRate'
        }
        dedup_result_series = dedup_result_series.loc[rename_dict.keys()].rename(rename_dict)

        dedup_result_series[f'FinalHisat3nReads'] = int(dedup_result_series[f'FilteredReads']) - int(
            dedup_result_series[f'DuplicatedReads']) + \
                                                    (int(dedup_result_series[f'FilteredPairs']) - int(
                                                        dedup_result_series[f'DuplicatedPairs'])) * 2
        dedup_result_series.name = cell_id
    except pd.errors.EmptyDataError:
        # if a BAM file is empty, picard matrix is also empty
        dedup_result_series = pd.Series({f'FilteredReads': 0,
                                         f'DuplicatedReads': 0,
                                         f'FilteredPairs': 0,
                                         f'DuplicatedPairs': 0,
                                         f'FinalHisat3nReads': 0,
                                         f'DuplicationRate': 0}, name=cell_id)
    return dedup_result_series
