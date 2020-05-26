import pathlib
from collections import defaultdict

import pandas as pd
import pysam

METHYLATED_CHAR = 'H'
UNMETHYLATED_CHAR = 'h'


def read_mc_level(bismark_tag):
    m_c = bismark_tag.count(METHYLATED_CHAR)
    normal_c = bismark_tag.count(UNMETHYLATED_CHAR)
    total_c = m_c + normal_c
    if total_c == 0:
        return 0, 0
    else:
        read_mc_rate = m_c / total_c
        return read_mc_rate, total_c


def select_dna_reads(input_bam,
                     output_bam,
                     mc_rate_max_threshold=0.5,
                     cov_min_threshold=3):
    read_profile_dict = defaultdict(int)
    # init dict to make sure the series has something
    read_profile_dict[(50, 50)] = 0
    with pysam.AlignmentFile(input_bam) as f:
        with pysam.AlignmentFile(output_bam, header=f.header, mode='wb') as out_f:
            for read in f:
                bismark_tag = read.get_tag('XM')
                mc_rate, cov = read_mc_level(bismark_tag)
                read_profile_dict[(int(100 * mc_rate), cov)] += 1

                # split reads
                if (mc_rate > mc_rate_max_threshold) or (cov < cov_min_threshold):
                    continue
                out_f.write(read)
    with open(str(output_bam) + '.reads_profile.csv', 'w') as stat_f:
        stat_f.write('mc_rate,cov,count\n')
        for (mc_rate, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_rate},{cov},{count}\n')
    return


def summarize_select_dna_reads(output_dir,
                               mc_rate_max_threshold=0.5,
                               cov_min_threshold=3):
    bam_dir = pathlib.Path(output_dir) / 'dna_bam'
    all_profile_path = bam_dir / 'select_dna_reads.all_profile.csv'
    final_stat_path = bam_dir / 'select_dna_reads.stats.csv'

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
    total_stats_df.to_csv(all_profile_path, index=None)

    selected_reads = total_stats_df[
        (total_stats_df['cov'] >= cov_min_threshold)
        & (total_stats_df['mc_rate'] < mc_rate_max_threshold)]

    selected_reads = selected_reads.groupby('cell_id')['count'].sum()
    selected_ratio = selected_reads / total_stats_df.groupby('cell_id')['count'].sum()
    final_stat = pd.DataFrame({'FinalDNAReads': selected_reads, 'FinalDNAReadsRatio': selected_ratio})
    final_stat.index.name = 'cell_id'
    final_stat.to_csv(final_stat_path, header=True)
    return


def dna_reads_selection(output_dir,
                        mc_rate_max_threshold=0.5,
                        cov_min_threshold=3):
    output_dir = pathlib.Path(output_dir).absolute()
    bam_batch = pd.read_csv(output_dir / 'snakemake/bismark_bam_list.txt',
                            header=None, index_col=0, squeeze=True)
    with open(output_dir / 'snakemake/bismark_bam_list.txt', 'w') as f:
        # clear bam list, replace it with dna_reads bam paths
        pass
    bam_dir = output_dir / 'dna_bam'
    bam_dir.mkdir(exist_ok=True)

    for batch_id, sub_series in bam_batch.groupby(bam_batch):
        bam_dict = {pathlib.Path(i).name.split('.')[0]: i
                    for i in sub_series.index}
        total_rules = ''
        dna_bam_paths = []
        for i, (cell_id, bam_path) in enumerate(bam_dict.items()):
            uid = '-'.join(cell_id.split('-')[:-1])
            this_bam_dir = bam_dir / uid
            this_bam_dir.mkdir(exist_ok=True)
            dna_bam_path = this_bam_dir / f'{cell_id}.dna_reads.bam'
            stat_path = this_bam_dir / f'{cell_id}.dna_reads.bam.reads_profile.csv'
            rule_template = f"""
rule select_dna_{i}:
    input:
        "{bam_path}"
    output:
        "{dna_bam_path}"
    log:
        "{stat_path}"
    shell:
        'yap-internal select-dna-reads --input_bam {{input}} '
        '--output_bam {{output}} --mc_rate_max_threshold {mc_rate_max_threshold} '
        '--cov_min_threshold {cov_min_threshold}'

"""
            total_rules += rule_template
            dna_bam_paths.append(str(dna_bam_path))

        with open(output_dir / f'snakemake/snakefile_select_dna_{batch_id}', 'w') as f:
            f.write(f"""
rule dna:
    input:
        {dna_bam_paths}

{total_rules}
""")
        # update bam list so
        with open(output_dir / 'snakemake/bismark_bam_list.txt', 'a') as f:
            for dna_bam_path in dna_bam_paths:
                f.write(f'{dna_bam_path},{batch_id}\n')
    return
