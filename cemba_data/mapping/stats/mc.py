import pathlib

import pandas as pd

from .utilities import parse_trim_fastq_stats, parse_trim_fastq_stats_mct, \
    parse_bismark_report, parse_deduplicate_stat, \
    generate_allc_stats


def mc_mapping_stats(output_dir, config):
    """this may apply to single UID dir, so config is provided as parameter"""
    output_dir = pathlib.Path(output_dir).absolute()
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'bam'
    allc_dir = output_dir / 'allc'
    cell_stats = []
    cell_ids = [path.name.split('.')[0]
                for path in allc_dir.glob(f'*.allc.tsv.gz')]

    for cell_id in cell_ids:
        print(f'Parsing stats of {cell_id}.')
        total_stats = []
        for read_type in ['R1', 'R2']:
            if config['mode'] == 'mct':
                total_stats.append(
                    parse_trim_fastq_stats_mct(
                        fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.txt'))
            else:
                total_stats.append(
                    parse_trim_fastq_stats(
                        fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.tsv'))
            total_stats.append(
                parse_bismark_report(
                    bam_dir / f'{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt'))
            total_stats.append(
                parse_deduplicate_stat(
                    bam_dir / f'{cell_id}-{read_type}.trimmed_bismark_bt2.deduped.matrix.txt'
                ))
        cell_stats.append(pd.concat(total_stats))
    mapping_df = pd.DataFrame(cell_stats)
    mapping_df.index.name = 'cell_id'

    # add allc stats
    allc_df = generate_allc_stats(output_dir, config)
    final_df = pd.concat([mapping_df, allc_df], sort=True, axis=1)
    return final_df


def mc_additional_cols(final_df):
    """Additional columns for mC mapping summary"""
    final_df = final_df.copy()
    final_df['CellInputReadPairs'] = final_df['R1InputReads'].astype(int)  # == final_df['R2InputReads']
    if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
        cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
                                        for _, i in final_df.groupby('PCRIndex')])
        final_df['CellBarcodeRatio'] = cell_barcode_ratio

    final_df['FinalmCReads'] = final_df['R1FinalBismarkReads'] + final_df['R2FinalBismarkReads']
    return final_df
