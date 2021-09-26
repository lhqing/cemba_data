import pathlib

import pandas as pd
import pysam

from .mc import mc_mapping_stats
from .mct import _count_reads_by_rg_in_star_bam, \
    summary_rna_mapping, \
    summarize_select_dna_reads, \
    aggregate_feature_counts
from .m3c import m3c_mapping_stats


def _4m_mapping_stats(output_dir, config):
    """this may apply to single UID dir, so config is provided as parameter"""
    m3c_stats_df = m3c_mapping_stats(output_dir, config)
    select_dna_stats_df = summarize_select_dna_reads(output_dir, config)
    rna_stats_df = summary_rna_mapping(output_dir)
    final_df = pd.concat([m3c_stats_df, select_dna_stats_df, rna_stats_df], axis=1)
    return final_df


def _4m_additional_cols(final_df, output_dir):
    final_df = final_df.copy()
    final_df['CellInputReadPairs'] = final_df['R1InputReads'].astype(int)  # == final_df['R2InputReads']
    if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
        cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
                                        for _, i in final_df.groupby('PCRIndex')])
        final_df['CellBarcodeRatio'] = cell_barcode_ratio

    stats = pd.read_hdf(output_dir / 'TotalRNAData.h5', key='stats')
    final_df['GenesDetected'] = stats['GenesDetected']

    # calculate some mCT specific ratios
    final_df['DNAReadsYield'] = final_df['FinalDNAReads'] / (
            final_df['CellInputReadPairs'] * 2)
    final_df['RNAReadsYield'] = final_df['FinalRNAReads'] / final_df[
        'CellInputReadPairs']
    final_df['RNA/(DNA+RNA)'] = final_df['FinalRNAReads'].fillna(0) / (
            final_df['R1FinalBismarkReads'].fillna(0) + 1)
    return final_df
