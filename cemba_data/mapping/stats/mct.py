import pathlib

import pandas as pd
import pysam

from .mc import mc_mapping_stats


def _count_reads_by_rg_in_star_bam(bam_path):
    try:
        bam = pysam.AlignmentFile(bam_path)
    except ValueError:
        # empty bam file
        return

    cell_read_counts = {cell['ID']: 0 for cell in bam.header['RG']}

    for read in bam:
        cell = read.get_tag('RG')
        cell_read_counts[cell] += 1
    read_count = pd.Series(cell_read_counts, name='Reads')
    read_count.index.name = 'cell_id'
    return read_count


def summary_rna_mapping(output_dir):
    output_dir = pathlib.Path(output_dir)

    # summarize read counts for each cell before filter by mC rate
    total_star_mapped_reads = _count_reads_by_rg_in_star_bam(output_dir / 'rna_bam/TotalRNAAligned.filtered.bam')

    # feature count summary
    total_counts = pd.read_csv(output_dir / 'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary',
                               sep='\t', index_col=0).T
    total_counts.index = total_counts.index.map(lambda i: i.split(':')[-1])
    feature_count_summary = total_counts[['Assigned']].copy()
    feature_count_summary['FinalRNAReads'] = total_counts.sum(axis=1)
    feature_count_summary.columns = ['FinalCountedReads', 'FinalRNAReads']

    total_rna_stat = feature_count_summary.copy()
    total_rna_stat['RNAUniqueMapped'] = total_star_mapped_reads
    total_rna_stat['SelectedRNAReadsRatio'] = total_rna_stat['FinalRNAReads'] / total_rna_stat['RNAUniqueMapped']
    total_rna_stat.index.name = 'cell_id'
    return total_rna_stat


def summarize_select_dna_reads(output_dir,
                               config):
    bam_dir = pathlib.Path(output_dir) / 'dna_bam'
    all_profile_path = bam_dir / 'select_dna_reads.all_profile.csv'
    final_stat_path = bam_dir / 'select_dna_reads.stats.csv'

    mc_rate_max_threshold = config['mc_rate_max_threshold']
    cov_min_threshold = config['dna_cov_min_threshold']

    records = []
    select_dna_reads_stat_list = list(bam_dir.glob('*/*.reads_profile.csv'))
    for path in select_dna_reads_stat_list:
        try:
            _df = pd.read_csv(path)
        except pd.errors.EmptyDataError:
            # means the bam file is empty
            continue

        cell_id = path.name.split('.')[0]
        _df['cell_id'] = cell_id
        _df['mc_rate_max_threshold'] = mc_rate_max_threshold
        _df['cov_min_threshold'] = cov_min_threshold
        records.append(_df)
    total_stats_df = pd.concat(records)

    selected_reads = total_stats_df[
        (total_stats_df['cov'] >= cov_min_threshold)
        & (total_stats_df['mc_rate'] < mc_rate_max_threshold)]

    selected_reads = selected_reads.groupby('cell_id')['count'].sum()
    selected_ratio = selected_reads / total_stats_df.groupby('cell_id')['count'].sum()
    final_stat = pd.DataFrame({'FinalDNAReads': selected_reads, 'SelectedDNAReadsRatio': selected_ratio})
    final_stat.index.name = 'cell_id'
    return final_stat


def mct_mapping_stats(output_dir, config):
    """this may apply to single UID dir, so config is provided as parameter"""
    mc_stats_df = mc_mapping_stats(output_dir, config)
    select_dna_stats = summarize_select_dna_reads(output_dir, config)
    rna_stats_df = summary_rna_mapping(output_dir)
    final_df = pd.concat([mc_stats_df, select_dna_stats, rna_stats_df])
    return final_df
