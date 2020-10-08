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
    total_rna_stat['RNAUniqueMappedReads'] = total_star_mapped_reads
    total_rna_stat['SelectedRNAReadsRatio'] = total_rna_stat['FinalRNAReads'] / total_rna_stat['RNAUniqueMappedReads']
    total_rna_stat.index.name = 'cell_id'
    return total_rna_stat


def summarize_select_dna_reads(output_dir,
                               config):
    bam_dir = pathlib.Path(output_dir) / 'bam'
    mc_rate_max_threshold = float(config['mc_rate_max_threshold'])
    cov_min_threshold = float(config['dna_cov_min_threshold'])

    records = []
    select_dna_reads_stat_list = list(bam_dir.glob('*.reads_profile.csv'))
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
    select_dna_stats_df = summarize_select_dna_reads(output_dir, config)
    rna_stats_df = summary_rna_mapping(output_dir)
    final_df = pd.concat([mc_stats_df, select_dna_stats_df, rna_stats_df], axis=1)
    return final_df


def aggregate_feature_counts(output_dir):
    output_dir = pathlib.Path(output_dir)
    cell_data = []

    count_paths = list(output_dir.glob('*/rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv'))
    if len(count_paths) == 0:
        return

    data = None
    for path in count_paths:
        data = pd.read_csv(path, sep='\t', index_col=0, comment='#')
        cell_data.append(data.iloc[:, 5:])
    cell_data = pd.concat(cell_data, axis=1, sort=True)
    cell_data.columns = cell_data.columns.str.split(':').str[1]

    # all count table should have the same info, so read the last one
    # chr, start, end, strand, length
    gene_info = data.iloc[:, :5]
    with pd.HDFStore(output_dir / 'TotalRNAData.h5', mode='w', complevel=5) as hdf:
        hdf['data'] = cell_data.T  # cell by gene
        hdf['gene'] = gene_info
        hdf['stats'] = pd.DataFrame({'GenesDetected': (cell_data > 0).sum()})
    return


def mct_additional_cols(final_df, output_dir):
    final_df = final_df.copy()
    final_df['CellInputReadPairs'] = final_df['R1InputReads'].astype(int)  # == final_df['R2InputReads']
    if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
        cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
                                        for _, i in final_df.groupby('PCRIndex')])
        final_df['CellBarcodeRatio'] = cell_barcode_ratio

    stats = pd.read_hdf(output_dir / 'TotalRNAData.h5', key='stats')
    final_df['GenesDetected'] = stats['GenesDetected']
    return final_df
