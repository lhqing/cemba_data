import pathlib
import subprocess

import pandas as pd

from .utilities import parse_trim_fastq_stats, generate_allc_stats


def _m3c_bam_unique_read_counts(bam_path):
    command = f"samtools view {bam_path} | awk -F ':N:0' '{{print $1}}' | sort | wc -l"
    p = subprocess.run(command, shell=True, check=True, encoding='utf8', stdout=subprocess.PIPE)
    unique_reads = p.stdout.strip()
    return unique_reads


def _m3c_count_bams(bam_dir, cell_id, read_type):
    bam_path_dict = {
        f'{read_type}UniqueMapped': bam_dir / f'{cell_id}-{read_type}.two_mapping.filter.bam',
        f'{read_type}Dedupped': bam_dir / f'{cell_id}-{read_type}.two_mapping.deduped.bam',
    }
    read_counts = {name: _m3c_bam_unique_read_counts(path)
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
