import logging
import pathlib
import shutil

import pandas as pd

import cemba_data

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def bismark_mapping(
        output_dir,
        bismark_reference,
        read_min=0,
        read_max=1e7,
        r1_adapter='AGATCGGAAGAGCACACGTCTGAAC',
        r1_left_cut=10,
        r1_right_cut=10,
        r2_adapter='AGATCGGAAGAGCGTCGTGTAGGGA',
        r2_left_cut=10,
        r2_right_cut=10,
        split_batch=True,
        batch_size=32,
        demultiplex_result=None,
        unmapped_fastq=False,
        batch_id=None,
        testing=False):
    # setup directories
    output_dir = pathlib.Path(output_dir).absolute()
    fastq_dir = output_dir / 'fastq'
    if not fastq_dir.exists():
        raise FileNotFoundError('Output dir do not have fastq/ sub dir.')
    snakemake_dir = output_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True)
    bam_dir = output_dir / 'bam'
    bam_dir.mkdir(exist_ok=True)
    read_min = max(0, read_min)

    unmapped_param_str = '--un' if unmapped_fastq else ''

    if demultiplex_result is None:
        demultiplex_result = pd.read_csv(fastq_dir / 'demultiplex.stats.csv',
                                         index_col='cell_id')

        # make sure snakemake_dir is empty
        if snakemake_dir.exists():
            shutil.rmtree(snakemake_dir)
            snakemake_dir.mkdir(exist_ok=True)

        if testing:
            if not isinstance(testing, int):
                testing = 10
            # testing X cell only
            demultiplex_result = demultiplex_result.iloc[:min(testing, demultiplex_result.shape[0]), :]

        if split_batch and demultiplex_result.shape[0] > batch_size:
            cell_chunks = (demultiplex_result[i:i + batch_size]
                           for i in range(0, demultiplex_result.shape[0], batch_size))
            for i, chunk in enumerate(cell_chunks):
                bismark_mapping(output_dir=output_dir,
                                bismark_reference=bismark_reference,
                                read_min=read_min,
                                read_max=read_max,
                                r1_adapter=r1_adapter,
                                r1_left_cut=r1_left_cut,
                                r1_right_cut=r1_right_cut,
                                r2_adapter=r2_adapter,
                                r2_left_cut=r2_left_cut,
                                r2_right_cut=r2_right_cut,
                                split_batch=split_batch,
                                batch_size=batch_size,
                                demultiplex_result=chunk,
                                unmapped_fastq=unmapped_fastq,
                                batch_id=i,
                                testing=testing)
            return
        else:
            batch_id = 0

    print('Input cells:', demultiplex_result.shape[0])
    demultiplex_result = demultiplex_result.query(f'{read_min} < CellInputReadPairs < {read_max}').copy()
    print('After CellInputReadPairs filter:', demultiplex_result.shape[0])
    total_rule = ''
    final_bams = []

    for i, (cell_id, row) in enumerate(demultiplex_result.iterrows()):
        cell_dir = fastq_dir / row['UID']
        this_dir = bam_dir / row['UID']
        this_dir.mkdir(exist_ok=True)
        picard_temp = this_dir / 'temp'
        picard_temp.mkdir(exist_ok=True)

        # R1
        r1_raw_fq = cell_dir / f'{cell_id}-R1.fq.gz'
        r1_trimmed_fq = this_dir / f'{cell_id}-R1.trimmed.fq.gz'
        r1_trim_stats = this_dir / f'{cell_id}-R1.trim_stats.tsv'
        # output name based on bismark behavior
        r1_raw_bam = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2.bam'
        r1_bismark_stats = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2_SE_report.txt'
        r1_filtered_bam = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2.filtered.bam'
        r1_sorted_bam = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2.sorted.bam'
        r1_dedup_bam = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2.final.bam'
        r1_dedup_stats = this_dir / f'{cell_id}-R1.trimmed_bismark_bt2.dedup.matrix.txt'

        # R2
        r2_raw_fq = cell_dir / f'{cell_id}-R2.fq.gz'
        r2_trimmed_fq = this_dir / f'{cell_id}-R2.trimmed.fq.gz'
        r2_trim_stats = this_dir / f'{cell_id}-R2.trim_stats.tsv'
        # output name based on bismark behavior
        r2_raw_bam = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2.bam'
        r2_bismark_stats = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2_SE_report.txt'
        r2_filtered_bam = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2.filtered.bam'
        r2_sorted_bam = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2.sorted.bam'
        r2_dedup_bam = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2.final.bam'
        r2_dedup_stats = this_dir / f'{cell_id}-R2.trimmed_bismark_bt2.dedup.matrix.txt'

        # final bam
        final_bam = this_dir / f'{cell_id}.final.bam'

        # rule for bismark mapping
        rule_template = f"""
rule trim_r1_{i}:
    input:
        "{r1_raw_fq}"
    output:
        "{r1_trimmed_fq}"
    log:
        "{r1_trim_stats}"
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r1_adapter} {{input}} 2> {{log}} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r1_left_cut} -u -{r1_right_cut} -m 30 "
        "-o {{output}} - >> {{log}}"

rule bismark_r1_{i}:
    input:
        "{r1_trimmed_fq}"
    output:
        temp("{r1_raw_bam}")
    threads: 
        3
    resources:
        mem_mb=16000
    log:
        "{r1_bismark_stats}"
    shell:
        "bismark {bismark_reference} {unmapped_param_str} --bowtie2 {{input}} --pbat -o {this_dir} --temp_dir {this_dir}"

rule filter_r1_bam_{i}:
    input:
        "{r1_raw_bam}"
    output:
        temp("{r1_filtered_bam}")
    shell:
        "samtools view -b -h -q 10 -o {{output}} {{input}}"

rule sort_r1_bam_{i}:
    input:
        "{r1_filtered_bam}"
    output:
        temp("{r1_sorted_bam}")
    shell:
        "samtools sort -o {{output}} {{input}}"

rule dedup_r1_bam_{i}:
    input:
        "{r1_sorted_bam}"
    output:
        temp("{r1_dedup_bam}")
    log:
        "{r1_dedup_stats}"
    shell:
        "picard MarkDuplicates I={{input}} O={{output}} M={{log}} REMOVE_DUPLICATES=true TMP_DIR={picard_temp}"

rule trim_r2_{i}:
    input:
        "{r2_raw_fq}"
    output:
        temp("{r2_trimmed_fq}")
    log:
        "{r2_trim_stats}"
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r2_adapter} {{input}} 2> {{log}} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r2_left_cut} -u -{r2_right_cut} -m 30 "
        "-o {{output}} - >> {{log}}"

rule bismark_r2_{i}:
    input:
        "{r2_trimmed_fq}"
    output:
        temp("{r2_raw_bam}")
    threads: 
        3
    resources:
        mem_mb=16000
    log:
        "{r2_bismark_stats}"
    shell:
        "bismark {bismark_reference} {unmapped_param_str} --bowtie2 {{input}} -o {this_dir} --temp_dir {this_dir}"

rule filter_r2_bam_{i}:
    input:
        "{r2_raw_bam}"
    output:
        temp("{r2_filtered_bam}")
    shell:
        "samtools view -b -h -q 10 -o {{output}} {{input}}"

rule sort_r2_bam_{i}:
    input:
        "{r2_filtered_bam}"
    output:
        temp("{r2_sorted_bam}")
    shell:
        "samtools sort -o {{output}} {{input}}"

rule dedup_r2_bam_{i}:
    input:
        "{r2_sorted_bam}"
    output:
        temp("{r2_dedup_bam}")
    log:
        "{r2_dedup_stats}"
    shell:
        "picard MarkDuplicates I={{input}} O={{output}} M={{log}} REMOVE_DUPLICATES=true TMP_DIR={picard_temp}"

rule merge_bam_{i}:
    input:
        ["{r1_dedup_bam}", "{r2_dedup_bam}"]
    output:
        "{final_bam}"
    shell:
        "samtools merge -f {{output}} {{input}}"

"""

        total_rule += rule_template
        final_bams.append(str(final_bam))

    with open(snakemake_dir / f'snakefile_bismark_mapping_{batch_id}', 'w') as f:
        f.write(f"""
rule final:
    input:
        {final_bams}
    output:
        touch("{output_dir}/snakemake/bismark_mapping_done_{batch_id}")

{total_rule}
""")
    # dump BAM paths into a file
    with open(snakemake_dir / 'bismark_bam_list.txt', 'a') as f:
        for bam_path in final_bams:
            f.write(f'{bam_path},{batch_id}\n')
    return


def _parse_trim_fastq_stats(stat_path):
    # example trim fastq stats
    """
status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
0	OK	1490	213724	0	0	0	1490	4	0	213712
1	status	in_reads	in_bp	too_short	too_long	too_many_n	out_reads	w/adapters	qualtrim_bp	out_bp
2	OK	1490	213712	0	0	0	1482	0	1300	182546
"""
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    trim_stats = pd.read_csv(stat_path, sep='\t')
    trim_stats = trim_stats.iloc[[0, 2], :].reset_index()  # skip the duplicated title row

    data = pd.Series({
        f'{read_type}InputReads': trim_stats['in_reads'][0],
        f'{read_type}InputReadsBP': trim_stats['in_bp'][0],
        f'{read_type}WithAdapters': trim_stats['w/adapters'][0],
        f'{read_type}QualTrimBP': trim_stats['qualtrim_bp'][1],
        f'{read_type}TrimmedReads': trim_stats['out_reads'][1],
        f'{read_type}TrimmedReadsBP': trim_stats['out_bp'][1],
        f'{read_type}TrimmedReadsRate': int(trim_stats['out_reads'][1]) / int(trim_stats['in_reads'][0])
    }, name=cell_id)
    return data


def _parse_bismark_report(stat_path):
    """
    parse bismark output
    """
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    term_dict = {
        'Number of alignments with a unique best hit from the different alignments': f'{read_type}UniqueMappedReads',
        'Mapping efficiency': f'{read_type}MappingRate',
        'Sequences with no alignments under any condition': f'{read_type}UnmappedReads',
        'Sequences did not map uniquely': f'{read_type}UnuniqueMappedReads',
        'CT/CT': f'{read_type}OT',
        'CT/GA': f'{read_type}OB',
        'GA/CT': f'{read_type}CTOT',
        'GA/GA': f'{read_type}CTOB',
        "Total number of C's analysed": f'{read_type}TotalC',
        'C methylated in CpG context': f'{read_type}TotalmCGRate',
        'C methylated in CHG context': f'{read_type}TotalmCHGRate',
        'C methylated in CHH context': f'{read_type}TotalmCHHRate'}

    with open(stat_path) as rep:
        report_dict = {}
        for line in rep:
            try:
                start, rest = line.split(':')
            except ValueError:
                continue  # more or less than 2 after split
            try:
                report_dict[term_dict[start]] = rest.strip().split('\t')[0].strip('%')
            except KeyError:
                pass
    return pd.Series(report_dict, name=cell_id)


def _parse_deduplicate_stat(stat_path):
    *cell_id, read_type = pathlib.Path(stat_path).name.split('.')[0].split('-')
    cell_id = '-'.join(cell_id)
    try:
        dedup_result_series = pd.read_csv(stat_path, comment='#', sep='\t').T[0]
        rename_dict = {
            'UNPAIRED_READS_EXAMINED': f'{read_type}MAPQFilteredReads',
            'UNPAIRED_READ_DUPLICATES': f'{read_type}DuplicatedReads',
            'PERCENT_DUPLICATION': f'{read_type}DuplicationRate'
        }
        dedup_result_series = dedup_result_series.loc[rename_dict.keys()].rename(rename_dict)

        dedup_result_series[f'{read_type}FinalBismarkReads'] = dedup_result_series[f'{read_type}MAPQFilteredReads'] - \
                                                               dedup_result_series[f'{read_type}DuplicatedReads']
        dedup_result_series.name = cell_id
    except pd.errors.EmptyDataError:
        # if a BAM file is empty, picard matrix is also empty
        dedup_result_series = pd.Series({f'{read_type}MAPQFilteredReads': 0,
                                         f'{read_type}DuplicatedReads': 0,
                                         f'{read_type}FinalBismarkReads': 0,
                                         f'{read_type}DuplicationRate': 0}, name=cell_id)
    return dedup_result_series


def bismark_mapping_stats(output_dir):
    output_dir = pathlib.Path(output_dir).absolute()
    bam_dir = output_dir / 'bam'
    cell_stats = []
    cell_ids = [path.name.split('.')[0]
                for path in bam_dir.glob('*/*.final.bam')]

    for cell_id in cell_ids:
        print(f'Parsing stats of {cell_id}.')
        uid = '-'.join(cell_id.split('-')[:-1])
        total_stats = []
        for read_type in ['R1', 'R2']:
            total_stats.append(
                _parse_trim_fastq_stats(
                    bam_dir / f'{uid}/{cell_id}-{read_type}.trim_stats.tsv'))
            total_stats.append(
                _parse_bismark_report(
                    bam_dir / f'{uid}/{cell_id}-{read_type}.trimmed_bismark_bt2_SE_report.txt'))
            total_stats.append(
                _parse_deduplicate_stat(
                    bam_dir / f'{uid}/{cell_id}-{read_type}.trimmed_bismark_bt2.dedup.matrix.txt'
                ))
        cell_stats.append(pd.concat(total_stats))
    final_df = pd.DataFrame(cell_stats)
    final_df.index.name = 'cell_id'
    final_df.to_csv(bam_dir / 'bismark_mapping_stats.csv')
    return
