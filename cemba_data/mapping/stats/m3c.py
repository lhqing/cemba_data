import pathlib

import pandas as pd
from pysam import AlignmentFile

from .utilities import parse_trim_fastq_stats, generate_allc_stats


def _m3c_bam_unique_read_counts(bam_path, read_type_int):
    unique_reads = set()
    with AlignmentFile(bam_path) as bam:
        for read in bam:
            unique_reads.add(read.query_name.split(f'_{read_type_int}:N:0:')[0])
    return len(unique_reads)


def _m3c_count_bams(bam_dir, cell_id, read_type):
    bam_path_dict = {
        f'{read_type}UniqueMappedReads': bam_dir / f'{cell_id}-{read_type}.two_mapping.filter.bam',
        f'{read_type}DeduppedReads': bam_dir / f'{cell_id}-{read_type}.two_mapping.deduped.bam',
    }
    read_counts = {name: _m3c_bam_unique_read_counts(path, 1 if read_type == 'R1' else 2)
                   for name, path in bam_path_dict.items()}
    return pd.Series(read_counts, name=cell_id)


def m3c_mapping_stats(output_dir, config):
    """this may apply to single UID dir, so config is provided as parameter"""
    output_dir = pathlib.Path(output_dir).absolute()
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'bam'
    hic_dir = output_dir / 'hic'
    cell_stats = []
    cell_ids = [path.name.split('.')[0]
                for path in bam_dir.glob('*.3C.sorted.bam')]

    for cell_id in cell_ids:
        total_stats = []  # list of series
        for read_type in ['R1', 'R2']:
            # fastq reads
            total_stats.append(
                parse_trim_fastq_stats(
                    fastq_dir / f'{cell_id}-{read_type}.trimmed.stats.tsv'))
            # bam reads
            total_stats.append(
                _m3c_count_bams(bam_dir, cell_id, read_type)
            )
        # contacts
        contact_counts = pd.read_csv(hic_dir / f'{cell_id}.3C.contact.tsv.gz.counts.txt',
                                     header=None, index_col=0, squeeze=True)
        contact_counts.name = cell_id
        total_stats.append(contact_counts)

        cell_stats.append(pd.concat(total_stats))
    total_df = pd.DataFrame(cell_stats)

    # add allc stats
    allc_df = generate_allc_stats(output_dir, config)
    final_df = pd.concat([total_df, allc_df], sort=True, axis=1)
    return final_df


def m3c_additional_cols(final_df):
    final_df['FinalmCReads'] = final_df['R1DeduppedReads'] + final_df['R2DeduppedReads']
    final_df['CellInputReadPairs'] = final_df['R1InputReads']
    # use % to be consistent with others
    final_df['R1MappingRate'] = final_df['R1UniqueMappedReads'] / final_df['R1TrimmedReads'] * 100
    final_df['R2MappingRate'] = final_df['R2UniqueMappedReads'] / final_df['R2TrimmedReads'] * 100
    final_df['R1DuplicationRate'] = (1 - final_df['R1DeduppedReads'] / final_df['R1UniqueMappedReads']) * 100
    final_df['R2DuplicationRate'] = (1 - final_df['R2DeduppedReads'] / final_df['R2UniqueMappedReads']) * 100

    if 'PCRIndex' in final_df.columns:  # plate info might not exist if the cell name is abnormal
        cell_barcode_ratio = pd.concat([(i['CellInputReadPairs'] / i['CellInputReadPairs'].sum())
                                        for _, i in final_df.groupby('PCRIndex')])
        final_df['CellBarcodeRatio'] = cell_barcode_ratio

    final_df['TotalContacts'] = final_df[
        ['CisShortContact', 'CisLongContact', 'TransContact']].sum(axis=1)
    final_df['CisShortRatio'] = final_df['CisShortContact'] / final_df['TotalContacts']
    final_df['CisLongRatio'] = final_df['CisLongContact'] / final_df['TotalContacts']
    final_df['TransRatio'] = final_df['TransContact'] / final_df['TotalContacts']
    return final_df
