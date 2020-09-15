import pathlib

import pandas as pd

from .utilities import parse_trim_fastq_stats, parse_bismark_report, parse_deduplicate_stat, generate_allc_stats


def snmc_mapping_stats(output_dir, config):
    """this may apply to single UID dir, so config is provided as parameter"""
    output_dir = pathlib.Path(output_dir).absolute()
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'bam'
    cell_stats = []
    cell_ids = [path.name.split('.')[0]
                for path in bam_dir.glob('*.final.bam')]

    for cell_id in cell_ids:
        print(f'Parsing stats of {cell_id}.')
        total_stats = []
        for read_type in ['R1', 'R2']:
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
