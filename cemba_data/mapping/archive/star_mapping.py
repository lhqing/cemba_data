import logging
import pathlib

import pandas as pd
import pysam

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def star_mapping(output_dir, star_reference, gtf_path, threads=30,
                 mc_rate_min_threshold=0.9, cov_min_threshold=3, remove_raw_bam=True,
                 feature_type='gene', id_type='gene_id'):
    output_dir = pathlib.Path(output_dir).absolute()
    cell_batch = pd.read_csv(output_dir / 'snakemake/bismark_bam_list.txt',
                             header=None, index_col=0, squeeze=True)
    # each batch will map together
    cell_batch.index = cell_batch.index.map(lambda path: pathlib.Path(path).name.split('.')[0])
    rna_dir = output_dir / 'rna_bam'
    rna_dir.mkdir(exist_ok=True)

    for batch_id, sub_series in cell_batch.groupby(cell_batch):
        total_rule = ''
        output_files = []
        star_size = 25
        for rule_id, start in enumerate(range(0, len(sub_series), star_size)):
            this_series = sub_series[start:start + star_size]

            trimmed_fastqs = []
            samples = this_series.index
            for i, cell_id in enumerate(this_series.index):
                uid = '-'.join(cell_id.split('-')[:-1])
                trimmed_fastq = output_dir / f'bam/{uid}/{cell_id}-R1.trimmed.fq.gz'
                trimmed_fastqs.append(str(trimmed_fastq))
            this_output_dir = output_dir / f'rna_bam/{batch_id}'
            this_output_dir.mkdir(exist_ok=True)
            output_prefix = this_output_dir / f'Total_{rule_id}_'
            read_files_in_str = ','.join(trimmed_fastqs)
            samples_str = ' , ID:'.join(samples)

            filtered_bam = f'{output_prefix}filtered.bam'
            rna_bam = f'{output_prefix}selected_rna.bam'
            feature_count_table = f'{output_prefix}feature_count.tsv'
            if remove_raw_bam:
                raw_bam_str = f'temp("{output_prefix}Aligned.out.bam")'
            else:
                raw_bam_str = f'"{output_prefix}Aligned.out.bam"'

            rule = f"""
rule star_{batch_id}_{rule_id}:
    input:
        {trimmed_fastqs}
    output:
        {raw_bam_str}
    log:
        "{output_prefix}Log.final.out"
    threads:
        {threads}
    shell:
        'STAR --runThreadN {{threads}} '
        '--genomeDir {star_reference} '
        '--alignEndsType EndToEnd '
        '--genomeLoad NoSharedMemory '
        '--outSAMstrandField intronMotif '
        '--outSAMtype BAM Unsorted '
        '--outSAMunmapped None '
        '--outSAMattributes NH HI AS NM MD '
        '--sjdbOverhang 100 '
        '--outFilterType BySJout '
        '--outFilterMultimapNmax 20 '
        '--alignSJoverhangMin 8 '
        '--alignSJDBoverhangMin 1 '
        '--outFilterMismatchNmax 999 '
        '--outFilterMismatchNoverLmax 0.04 '
        '--alignIntronMin 20 '
        '--alignIntronMax 1000000 '
        '--alignMatesGapMax 1000000 '
        '--outFileNamePrefix {output_prefix} '
        '--readFilesIn {read_files_in_str} '
        '--readFilesCommand gzip -cd '
        '--outSAMattrRGline ID:{samples_str}'

rule filter_bam_{batch_id}_{rule_id}:
    input:
        "{output_prefix}Aligned.out.bam"
    output:
        "{filtered_bam}"
    threads:
        {threads}
    shell:
        "samtools sort -@ {{threads}} -m 2G {{input}} | samtools view -bh -q 10 -o {{output}} -"

rule select_rna_{batch_id}_{rule_id}:
    input:
        "{filtered_bam}"
    output:
        "{rna_bam}"
    shell:
        'yap-internal select-rna-reads ' \
        '--input_bam {{input}} ' \
        '--output_bam {{output}} ' \
        '--mc_rate_min_threshold {mc_rate_min_threshold} ' \
        '--cov_min_threshold {cov_min_threshold} '

rule feature_count_{batch_id}_{rule_id}:
    input:
        "{rna_bam}"
    output:
        "{feature_count_table}"
    threads:
        {threads}
    shell:
        'featureCounts -t {feature_type} -g {id_type} ' \
        '-a {gtf_path} -o {{output}} --byReadGroup -T {{threads}} {{input}}'
"""
            total_rule += rule
            output_files.append(feature_count_table)
        final = f"""
include: "{output_dir}/snakemake/snakefile_bismark_mapping_{batch_id}"

rule rna:
    input:
        {output_files}
    output:
        touch("{output_dir}/snakemake/star_mapping_done_{batch_id}")

{total_rule}
"""
        with open(output_dir / f'snakemake/snakefile_star_mapping_{batch_id}', 'w') as f:
            f.write(final)
    return


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
    read_count.to_csv(str(bam_path) + '.cell_stats.csv', header=True)
    return


def summary_rna_mapping(output_dir):
    output_dir = pathlib.Path(output_dir)

    # summarize read counts for each cell before filter by mC rate
    for path in output_dir.glob('rna_bam/*/Total*filtered.bam'):
        _count_reads_by_rg_in_star_bam(path)
    total_star_mapped_reads = pd.concat([pd.read_csv(path, index_col='cell_id')
                                         for path in output_dir.glob('rna_bam/*/Total*filtered.bam.cell_stats.csv')])
    total_star_mapped_reads.columns = ['RNAUniqueMapped']

    # feature count summary
    feature_count_summary_list = list(output_dir.glob('rna_bam/*/*_feature_count.tsv.summary'))
    total_counts = pd.concat([pd.read_csv(p, sep='\t', index_col=0).T for p in feature_count_summary_list])
    total_counts.index = total_counts.index.map(lambda i: i.split(':')[-1])
    feature_count_summary = total_counts[['Assigned']].copy()
    feature_count_summary['FinalRNAReads'] = total_counts.sum(axis=1)
    feature_count_summary.columns = ['FinalCountedReads', 'FinalRNAReads']

    total_rna_stat = pd.concat([total_star_mapped_reads, feature_count_summary], axis=1)
    total_rna_stat['SelectedRNAReadsRatio'] = total_rna_stat['FinalRNAReads'] / total_rna_stat['RNAUniqueMapped']
    total_rna_stat.index.name = 'cell_id'
    total_rna_stat.to_csv(output_dir / 'rna_bam/star_mapping_stats.csv')
    return
