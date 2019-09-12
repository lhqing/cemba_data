import collections
import logging
import os
import pathlib

import pandas as pd

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def print_default_configuration(out_path=None):
    """
    Print default .ini config file or save to out_path
    """
    with open(os.path.dirname(__file__) + '/mapping_config.ini') as f:
        configs = f.readlines()
    if out_path is not None:
        with open(out_path, 'w') as f:
            f.writelines(configs)
    else:
        for line in configs:
            print(line, end='')
    return


def summary_pipeline_stat(out_dir):
    """
    Combine all statistics and do some additional computation.

    Parameters
    ----------
    out_dir
        Pipeline universal directory

    Returns
    -------
    total_meta
        dataframe contain every parts of metadata.
        Each row is a cell, some metadata are reduced in some way, such as lanes.
    """
    out_dir = pathlib.Path(out_dir)

    # merge all stats dataframe
    stat_dict = collections.defaultdict(list)
    for f in out_dir.glob('**/stats/*tsv.gz'):
        df = pd.read_csv(f, sep='\t')
        stat_dict[f.name.split('.')[0]].append(df)
    # result df is a dict contain 5 concat dfs for raw metadata
    result_dfs = {}
    for k, v in stat_dict.items():
        result_dfs[k] = pd.concat(v, ignore_index=True, sort=True)
    # save all concatenated stats into output_dir level stats folder
    total_stats = out_dir / 'stats'
    total_stats.mkdir(exist_ok=True)
    for k, v in result_dfs.items():
        k_path = total_stats / f'{k}.tsv.gz'
        v.to_csv(str(k_path), sep='\t', compression='gzip')

    # demultiplex stat, from cutadapt demultiplex step
    demultiplex_result = result_dfs['demultiplex_records'].groupby(['uid', 'index_name']) \
        .sum()[['TotalPair', 'Trimmed']]
    demultiplex_result.rename(columns={'TotalPair': 'MultiplexReadsTotal',
                                       'Trimmed': 'IndexReadsTotal'},
                              inplace=True)
    demultiplex_result['MultiplexReadsTotal'] *= 2
    demultiplex_result['IndexReadsTotal'] *= 2
    demultiplex_result['IndexReadsRatio'] = demultiplex_result['IndexReadsTotal'] / demultiplex_result[
        'MultiplexReadsTotal'] * 100

    # fastq trim stat
    fastq_trim_result = result_dfs['fastq_trim_result'].groupby(['uid', 'index_name']).sum()[
        ['in_bp', 'out_bp', 'out_reads', 'qualtrim_bp', 'too_short', 'w/adapters']]
    fastq_trim_result.rename(columns={'in_bp': 'IndexBpTotal',
                                      'out_bp': 'IndexTrimedBpTotal',
                                      'out_reads': 'IndexTrimedReadsTotal',
                                      'qualtrim_bp': 'ReadsQualTrimBpTotal',
                                      'too_short': 'ReadsLengthFilterTotal',
                                      'w/adapters': 'ReadsWithAdapterTotal'}, inplace=True)
    fastq_trim_result['TrimedReadsAveLength'] = fastq_trim_result['IndexTrimedBpTotal'] / fastq_trim_result[
        'IndexTrimedReadsTotal']
    fastq_trim_result['IndexTrimedReadsRatio'] = fastq_trim_result['IndexTrimedReadsTotal'] / demultiplex_result[
        'IndexReadsTotal']

    # bismark_mapping stat
    bismark_r1 = result_dfs['bismark_result'][result_dfs['bismark_result']['read_type'].apply(lambda i: '1' in i)] \
        .set_index(['uid', 'index_name'])[['CTOB', 'CTOT', 'mapping_rate', 'total_c',
                                           'total_reads', 'unique_map', 'unmap', 'ununique_map']]
    bismark_r1.rename(columns={'mapping_rate': 'R1MappedRatio',  # this mapping rate is unique mapping rate
                               'total_c': 'R1TotalC',
                               'total_reads': 'R1TrimmedReads',
                               'unique_map': 'R1UniqueMappedReads',
                               'unmap': 'R1UnmappedReads',
                               'ununique_map': 'R1UnuniqueMappedReads'},
                      inplace=True)
    bismark_r2 = result_dfs['bismark_result'][result_dfs['bismark_result']['read_type'].apply(lambda i: '2' in i)] \
        .set_index(['uid', 'index_name'])[['OB', 'OT', 'mapping_rate', 'total_c',
                                           'total_reads', 'unique_map', 'unmap', 'ununique_map']]
    bismark_r2.rename(columns={'mapping_rate': 'R2MappedRatio',
                               'total_c': 'R2TotalC',
                               'total_reads': 'R2TrimmedReads',
                               'unique_map': 'R2UniqueMappedReads',
                               'unmap': 'R2UnmappedReads',
                               'ununique_map': 'R2UnuniqueMappedReads'},
                      inplace=True)
    bismark_result = pd.concat([bismark_r1, bismark_r2], sort=True, axis=1)
    bismark_result['TotalUniqueMappedReads'] = \
        bismark_result['R1UniqueMappedReads'] + bismark_result['R2UniqueMappedReads']
    bismark_result['TotalMappedRatio'] = bismark_result['TotalUniqueMappedReads'] / bismark_result[[
        'R1TrimmedReads', 'R2TrimmedReads']].sum(axis=1)

    bam_result = result_dfs['bam_process_result'].groupby(['uid', 'index_name']) \
        .sum()[['UNPAIRED_READ_DUPLICATES', 'out_reads']] \
        .rename(columns={'UNPAIRED_READ_DUPLICATES': 'DupReads',
                         'out_reads': 'DeduppedReads'})
    bam_result['DeduppedRatio'] = bam_result['DeduppedReads'] / bam_result.sum(axis=1)

    # ALLC stat
    cov_df = result_dfs['stat_allc_total_result'] \
        .set_index(['uid', 'index_name', 'index'])['cov'] \
        .unstack('index') \
        .fillna(0).astype(int)
    ccc_cov = cov_df['CCC']
    cov_df = cov_df.groupby(lambda i: i[:2], axis=1) \
        .sum()[['CA', 'CC', 'CG', 'CT']] \
        .rename(columns={c: c + '_Cov' for c in ['CA', 'CC', 'CG', 'CT']})
    cov_df['CH_Cov'] = cov_df[['CA_Cov', 'CC_Cov', 'CT_Cov']].sum(axis=1)
    mc_df = result_dfs['stat_allc_total_result'] \
        .set_index(['uid', 'index_name', 'index'])['mc'] \
        .unstack('index').fillna(0).astype(int)
    ccc_mc = mc_df['CCC']
    mc_df = mc_df.groupby(lambda i: i[:2], axis=1) \
        .sum()[['CA', 'CC', 'CG', 'CT']] \
        .rename(columns={c: c + '_Mc' for c in ['CA', 'CC', 'CG', 'CT']})
    mc_df['CH_Mc'] = mc_df[['CA_Mc', 'CC_Mc', 'CT_Mc']].sum(axis=1)
    # add mc rate
    mc_df['CH_Rate'] = mc_df['CH_Mc'] / cov_df['CH_Cov']
    mc_df['CG_Rate'] = mc_df['CG_Mc'] / cov_df['CG_Cov']
    # add CCC rate and use it to add adj mc rate
    mc_df['CCC_Mc'] = ccc_mc
    mc_df['CCC_Cov'] = ccc_cov
    mc_df['CCC_Rate'] = ccc_mc / ccc_cov
    mc_df['CH_RateAdj'] = (mc_df['CH_Rate'] - mc_df['CCC_Rate']) / (1 - mc_df['CCC_Rate'])
    mc_df['CG_RateAdj'] = (mc_df['CG_Rate'] - mc_df['CCC_Rate']) / (1 - mc_df['CCC_Rate'])

    # concat total meta
    total_meta = pd.concat([demultiplex_result, fastq_trim_result, bismark_result,
                            bam_result, mc_df, cov_df],
                           sort=True, axis=1).dropna()

    # file paths
    allc_dict = {}
    for f in out_dir.glob('**/allc*tsv.gz'):
        if 'stats' in f.parent.name:
            continue
        _, uid, index_name = (f.name.split('.')[0].split('_'))
        allc_dict[(uid, index_name)] = str(f.absolute().name)
    total_meta['AllcPath'] = pd.Series(allc_dict)
    return total_meta


def pipeline(input_fastq_pattern, output_dir, config_path):
    # test environment
    from .test_environment import testing_mapping_installation
    testing_mapping_installation()

    # pipeline_fastq


    # pipeline_mc

    # if mct: pipeline_rna

    # pipeline mcds

    return 0
