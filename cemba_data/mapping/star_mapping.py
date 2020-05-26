import logging
import pathlib
import shlex
import subprocess

import pandas as pd
import pysam

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def star_mapping(output_dir, star_reference, threads=30,
                 mc_rate_min_threshold=0.9, cov_min_threshold=3):
    output_dir = pathlib.Path(output_dir).absolute()
    cell_batch = pd.read_csv(output_dir / 'snakemake/bismark_bam_list.txt',
                             header=None, index_col=0, squeeze=True)
    # each batch will map together
    cell_batch.index = cell_batch.index.map(lambda path: pathlib.Path(path).name.split('.')[0])
    rna_dir = output_dir / 'rna_bam'
    rna_dir.mkdir(exist_ok=True)

    for batch_id, sub_series in cell_batch.groupby(cell_batch):
        trimmed_fastqs = []
        samples = sub_series.index
        for i, cell_id in enumerate(sub_series.index):
            uid = '-'.join(cell_id.split('-')[:-1])
            trimmed_fastq = output_dir / f'bam/{uid}/{cell_id}-R1.trimmed.fq.gz'
            trimmed_fastqs.append(str(trimmed_fastq))
        this_output_dir = output_dir / f'rna_bam/{batch_id}'
        this_output_dir.mkdir(exist_ok=True)
        output_prefix = this_output_dir / f'Total_'
        read_files_in_str = ','.join(trimmed_fastqs)
        samples_str = ' , ID:'.join(samples)

        filtered_bam = f'{output_prefix}filtered.bam'
        rna_bam = f'{output_prefix}selected_rna.bam'

        rule = f"""
rule rna:
    input:
        "{rna_bam}"

rule star_{batch_id}:
    input:
        {trimmed_fastqs}
    output:
        temp("{output_prefix}Aligned.out.bam")
    log:
        "{output_prefix}Log.final.out"
    threads:
        {threads}
    shell:
        'STAR --runThreadN {{threads}} '
        '--genomeDir {star_reference} '
        '--alignEndsType EndToEnd '
        '--genomeLoad LoadAndKeep '
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

rule filter_bam_{batch_id}:
    input:
        "{output_prefix}Aligned.out.bam"
    output:
        "{filtered_bam}"
    threads:
        {threads}
    shell:
        "samtools sort -@ {{threads}} -m 2G {{input}} | samtools view -bh -q 10 -o {{output}} -"

rule select_rna_{batch_id}:
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

"""
        with open(output_dir / f'snakemake/snakefile_star_mapping_{batch_id}', 'w') as f:
            f.write(rule)
    return


def _count_rna_reads_per_cell(bam_path, split=False):
    bam = pysam.AlignmentFile(bam_path)
    cell_read_counts = {cell['ID']: 0 for cell in bam.header['RG']}

    for read in bam:
        cell = read.get_tag('RG')
        cell_read_counts[cell] += 1
    read_count = pd.Series(cell_read_counts, name='Reads')
    read_count.index.name = 'cell_id'
    read_count.to_csv(str(bam_path) + '.cell_stats.csv', header=True)

    if split:
        output_dir = pathlib.Path(bam_path).parent
        subprocess.run(shlex.split(f'samtools split {bam_path} -f {output_dir}/%\!.rna_reads.%.'),
                       check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    return


def _post_star_mapping_process(bam_dir):
    if not pathlib.Path(f'{bam_dir}/Total_selected_rna.bam').exists():
        return
    _count_rna_reads_per_cell(f'{bam_dir}/Total_filtered.bam')
    _count_rna_reads_per_cell(f'{bam_dir}/Total_selected_rna.bam', split=True)

    file_names = [
        'Total_filtered.bam', 'Total_Log.final.out', 'Total_Log.out',
        'Total_Log.progress.out', 'Total_selected_rna.bam', 'Total_SJ.out.tab'
    ]

    subprocess.run(['rm', '-f'] + [f'{bam_dir}/{i}' for i in file_names], check=True)

    bam_list = pathlib.Path(bam_dir).absolute().glob('*.rna_reads.bam')
    bam_series = pd.Series({bam_path.name.split('.')[0]: str(bam_path) for bam_path in bam_list})
    bam_series.to_csv(f'{bam_dir}/cell_rna_reads_bam.csv', header=False)
    return


def summary_rna_mapping(output_dir):
    output_dir = pathlib.Path(output_dir)
    batch_sub_dirs = [p for p in output_dir.glob('rna_bam/*') if (p.is_dir())]
    for sub_dir in batch_sub_dirs:
        _post_star_mapping_process(sub_dir)

    total_star_mapped_reads = pd.concat([pd.read_csv(path, index_col='cell_id')
                                         for path in output_dir.glob('rna_bam/*/Total_filtered.bam.cell_stats.csv')])
    total_star_mapped_reads.columns = ['RNAUniqueMapped']

    total_rna_reads = pd.concat([pd.read_csv(path, index_col='cell_id')
                                 for path in output_dir.glob('rna_bam/*/Total_selected_rna.bam.cell_stats.csv')])
    total_rna_reads.columns = ['FinalRNAReads']

    total_rna_stat = pd.concat([total_star_mapped_reads, total_rna_reads], axis=1)

    total_rna_stat['FinalRNAReadsRatio'] = total_rna_stat['FinalRNAReads'] / total_rna_stat['RNAUniqueMapped']
    total_rna_stat.to_csv(output_dir / 'rna_bam/star_mapping_stats.csv')
    return



