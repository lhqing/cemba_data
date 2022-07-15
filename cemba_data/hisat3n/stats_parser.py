"""
Each cell_parser function takes a path from single cell and return a series named by cell_id

"""

import pathlib

import numpy as np
import pandas as pd

from .stats_col_names import COL_NAMES


def cell_parser_hisat_summary(stat_path):
    """
    parse hisat3n summary file
    """
    cell_id = pathlib.Path(stat_path).name.split('.')[0]
    term_dict = {
        'Total pairs': 'ReadPairsMappedInPE',
        'Aligned concordantly or discordantly 0 time': 'PEUnmappableReadPairs',
        'Aligned concordantly 1 time': 'PEUniqueMappedReadPairs',
        'Aligned concordantly >1 times': 'PEMultiMappedReadPairs',
        'Aligned discordantly 1 time': 'PEDiscordantlyUniqueMappedReadPairs',
        'Total unpaired reads': 'ReadsMappedInSE',
        'Aligned 0 time': 'SEUnmappableReads',
        'Aligned 1 time': 'SEUniqueMappedReads',
        'Aligned >1 times': 'SEMultiMappedReads',
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
                report_dict[term_dict[start]] = rest.strip().split(
                    ' ')[0].strip('%')
            except KeyError:
                pass
        for v in term_dict.values():
            if v not in report_dict:
                report_dict[v] = 0

        report_dict = pd.Series(report_dict).astype(int)
        total_reads = report_dict[f'ReadPairsMappedInPE'] * 2 + report_dict[
            'ReadsMappedInSE']
        unique_mapped_reads = report_dict[f'PEUniqueMappedReadPairs'] * 2 + \
                              report_dict[f'PEDiscordantlyUniqueMappedReadPairs'] * 2 + \
                              report_dict[f'SEUniqueMappedReads']
        report_dict['UniqueMappedReads'] = unique_mapped_reads
        report_dict[f'UniqueMappingRate'] = round(unique_mapped_reads /
                                                  (total_reads + 0.00001) * 100)
        multi_mapped_reads = report_dict[f'PEMultiMappedReadPairs'] * 2 + \
                             report_dict[f'SEMultiMappedReads']
        report_dict['MultiMappedReads'] = multi_mapped_reads
        report_dict[f'MultiMappingRate'] = round(multi_mapped_reads /
                                                 (total_reads + 0.00001) * 100)
        report_dict[f'OverallMappingRate'] = round(
            (unique_mapped_reads + multi_mapped_reads) / (total_reads + 0.00001) * 100)
    return pd.Series(report_dict, name=cell_id)


def cell_parser_picard_dedup_stat(stat_path):
    stat_path = pathlib.Path(stat_path)
    cell_id, *other_parts = stat_path.name.split('.')
    try:
        record = pd.read_csv(stat_path, comment='#', sep='\t').T[0]
        record['FinalReads'] = int(record['UNPAIRED_READS_EXAMINED']) - \
                               int(record['UNPAIRED_READ_DUPLICATES']) + \
                               (int(record['READ_PAIRS_EXAMINED']) - int(record['READ_PAIR_DUPLICATES'])) * 2
        record['DuplicatedReads'] = int(record['UNPAIRED_READ_DUPLICATES']) + \
                                    int(record['READ_PAIR_DUPLICATES']) * 2
        record['PCRDuplicationRate'] = record['FinalReads'] / \
                                       (record['FinalReads'] + record['DuplicatedReads'])
        record['PCRDuplicationRate'] = int(
            (1 - record['PCRDuplicationRate']) * 100)

        record.name = cell_id
    except pd.errors.EmptyDataError:
        # if a BAM file is empty, picard matrix is also empty
        record = pd.Series(
            {
                k: 0
                for p, k in COL_NAMES.keys()
                if p == 'parse_picard_dedup_stat'
            },
            dtype='int',
            name=cell_id)
    return record


def cell_parser_cutadapt_trim_stats(path):
    path = pathlib.Path(path)

    cell_id = path.name.split('.')[0]
    cell_records = pd.read_csv(path, sep='\t').T.squeeze()
    cell_records.name = cell_id
    return cell_records


def cell_parser_allc_count(path):
    path = pathlib.Path(path)
    cell_id = path.name.split('.')[0]

    allc_counts = pd.read_csv(path, index_col=0)

    if allc_counts.empty:
        return pd.Series([], dtype='O', name=cell_id)

    # remove contexts that contain "N"
    allc_counts = allc_counts.loc[allc_counts.index.map(
        lambda a: 'N' not in a)].copy()

    # find out which position is the C
    try:
        assert allc_counts.index.map(lambda a: len(a)).unique().size == 1
        c_pos = None
        for i in range(0, len(allc_counts.index[0])):
            if allc_counts.index.str[i].unique().size == 1:
                c_pos = i
        assert c_pos is not None
        assert c_pos != len(allc_counts.index[0])
    except AssertionError:
        raise AssertionError(f'Do not understand the mC context in {path}')

    # get mC context
    mc_context = pd.Series(allc_counts.index.str[c_pos + 1] == 'G').map({
        True: 'mCG', False: 'mCH'
    })
    if c_pos > 0:
        # NOMe
        nome_sites = pd.Series(allc_counts.index.str[c_pos - 1] == 'G').map({
            True: 'G', False: 'H'
        })
        mc_context = nome_sites + mc_context
    mc_context.index = allc_counts.index

    # generate cell records
    mc_context_sum = allc_counts.groupby(mc_context).sum()[['mc', 'cov']]

    is_ccc = allc_counts.index.str[c_pos:].map(
        lambda a: len(set(a)) == 1).values
    if c_pos > 0:
        is_ccc &= np.array(allc_counts.index.str[c_pos - 1] != 'G')

    try:
        ccc_mc, ccc_cov = np.ravel(allc_counts.loc[is_ccc, ['mc', 'cov']].values)
    except ValueError:
        ccc_mc, ccc_cov = 0, 0

    mc_context_sum = pd.concat([
        pd.DataFrame({
            'mCCC': {
                'mc': ccc_mc,
                'cov': ccc_cov
            }
        }).T, mc_context_sum
    ])
    mc_context_sum = mc_context_sum.astype('O')
    mc_context_sum['Frac'] = mc_context_sum['mc'] / (mc_context_sum['cov'] +
                                                     0.00001)
    mc_context_sum.rename(columns={'mc': 'mC', 'cov': 'Cov'}, inplace=True)

    cell_records = {}
    for (count_type, mc_type), count in mc_context_sum.unstack().items():
        cell_records[f'{mc_type}{count_type}'] = count
    cell_records = pd.Series(cell_records, name=cell_id, dtype='O')
    return cell_records


def cell_parser_reads_mc_frac_profile(path):
    # read parameters
    params = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                k, v = line.strip().split('=')
                params[k[1:]] = v
            else:
                break

    params['cov_min_threshold'] = int(params['cov_min_threshold'])
    # the mc_frac col in stats file is in %
    try:
        params['mc_rate_max_threshold'] = float(params['mc_rate_max_threshold'])
        if params['mc_rate_max_threshold'] <= 1:
            params['mc_rate_max_threshold'] *= 100
    except ValueError:
        pass
    try:
        params['mc_rate_min_threshold'] = float(params['mc_rate_min_threshold'])
        if params['mc_rate_min_threshold'] <= 1:
            params['mc_rate_min_threshold'] *= 100
    except ValueError:
        pass

    path = pathlib.Path(path)
    cell_id = path.name.split('.')[0]
    reads_profile = pd.read_csv(path, comment='#')
    mode = params['mode'].upper()
    if mode == 'DNA':
        selected_reads = reads_profile[
            (reads_profile['cov'] >= params['cov_min_threshold'])
            & (reads_profile['mc_frac'] <= params['mc_rate_max_threshold'])]
    else:
        selected_reads = reads_profile[
            (reads_profile['cov'] >= params['cov_min_threshold'])
            & (reads_profile['mc_frac'] >= params['mc_rate_min_threshold'])]

    selected_reads = selected_reads['count'].sum()
    selected_ratio = selected_reads / (reads_profile['count'].sum() + 0.0001)
    final_stat = pd.Series({f'Final{mode}Reads': selected_reads,
                            f'Selected{mode}ReadsRatio': selected_ratio},
                           dtype='O', name=cell_id)
    return final_stat


def cell_parser_feature_count_summary(path):
    result = pd.read_csv(path, sep='\t', index_col=0).squeeze()
    all_reads = result.sum()

    result.index.name = None
    result.name = result.name.split(':')[-1]  # cell id
    result['Unassigned_Total'] = result[result.index.str.startswith('Unassigned')].sum()
    result['AssignedRNAReadsRate'] = int(result['Assigned'] / (all_reads + 0.00001) * 100)
    return result


def cell_parser_call_chromatin_contacts(path):
    path = pathlib.Path(path)
    contact_stats = pd.read_csv(
        path,
        header=None,
        index_col=0).squeeze()
    contact_stats.name = path.name.split('.')[0]

    # add some calculated columns
    contact_stats['TotalCisContacts'] = contact_stats[contact_stats.index.str.startswith('cis')].sum()
    contact_stats['TotalTransContacts'] = contact_stats[contact_stats.index.str.startswith('trans')].sum()
    contact_stats['TotalMultiContacts'] = contact_stats[contact_stats.index.str.endswith('_multi')].sum()
    contact_stats['CisContactsRatio'] = contact_stats['TotalCisContacts'] / \
                                        (contact_stats['mapped_frag'] + 0.00001)
    contact_stats['TransContactsRatio'] = contact_stats['TotalTransContacts'] / \
                                          (contact_stats['mapped_frag'] + 0.00001)
    contact_stats['MultiContactsRatio'] = contact_stats['TotalMultiContacts'] / \
                                          (contact_stats['mapped_frag'] + 0.00001)
    return contact_stats


def parse_single_stats_set(path_pattern, parser, prefix=''):
    """
    Parse all the stats files in the path_pattern and return a dataframe
    with each cell id as index and the stats as columns.

    Parameters
    ----------
    path_pattern :
        Path pattern to the stats files.
    parser :
        A function that takes a path and returns a pandas series for one cell.
    prefix :
        Prefix to add to the column names.

    Returns
    -------
    pd.DataFrame
    """
    detail_stats_dir = pathlib.Path('detail_stats/')
    detail_stats_dir.mkdir(exist_ok=True)

    stats_paths = list(pathlib.Path().glob(path_pattern))

    records = []
    empty_index = []
    for path in stats_paths:
        try:
            record = parser(path)
            if record.size == 0:
                empty_index.append(record.name)
            else:
                records.append(record)
        except BaseException as e:
            print(f'Got error when reading {path} with {parser}')
            raise e
    stats_df = pd.DataFrame(records)
    if len(empty_index) > 0:
        use_index = pd.Index(stats_df.index.tolist() + empty_index)
        stats_df = stats_df.reindex(use_index)  # still record empty entries, but values will be nan

    # before rename, save raw stats into detail_stats/
    stats_df.to_csv(f'detail_stats/{prefix}.{parser.__name__}.csv')

    # rename columns
    rename_dict = {}
    rename_functions = set([k[0] for k in COL_NAMES.keys()])
    if parser.__name__ not in rename_functions:
        # do not rename function output if the function name not in COL_NAMES
        return stats_df

    for (parser_name, original_name), new_name in COL_NAMES.items():
        if parser_name == parser.__name__:
            if new_name == '':
                rename_dict[original_name] = original_name
            else:
                rename_dict[original_name] = new_name
    new_columns = stats_df.columns.map(rename_dict)
    if new_columns.isna().sum() != 0:
        print(
            f'These columns are ignored because they do not listed'
            f' in the COL_NAMES dict: {stats_df.columns[new_columns.isna()]}')
    # remove columns marked as DELETE in COL_NAMES
    stats_df.columns = new_columns
    stats_df = stats_df.iloc[:, new_columns != 'DELETE'].copy()

    # add prefix to stats_df columns
    stats_df.columns = stats_df.columns.map(lambda a: prefix + str(a))
    return stats_df
