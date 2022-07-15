from .stats_parser import *


def snmc_summary():
    """
    Generate snmC pipeline MappingSummary.csv.gz and save into cwd

    Returns
    -------
    pd.DataFrame
    """
    all_stats = []

    # fastq trimming stats
    df = parse_single_stats_set(f'fastq/*.trimmed.stats.txt',
                                cell_parser_cutadapt_trim_stats)
    all_stats.append(df)

    # hisat-3n mapping
    df = parse_single_stats_set(f'bam/*.hisat3n_dna_summary.txt',
                                cell_parser_hisat_summary)
    all_stats.append(df)

    # uniquely mapped reads dedup
    df = parse_single_stats_set(f'bam/*.unique_align.deduped.matrix.txt',
                                cell_parser_picard_dedup_stat, prefix='UniqueAlign')
    all_stats.append(df)

    # multi mapped reads dedup
    df = parse_single_stats_set(f'bam/*.multi_align.deduped.matrix.txt',
                                cell_parser_picard_dedup_stat, prefix='MultiAlign')
    all_stats.append(df)

    # allc count
    df = parse_single_stats_set(f'allc/*.allc.tsv.gz.count.csv',
                                cell_parser_allc_count)
    all_stats.append(df)

    # concatenate all stats
    all_stats = pd.concat(all_stats, axis=1)
    all_stats.index.name = 'cell'
    all_stats.to_csv(f'MappingSummary.csv.gz')
    return all_stats


def snmct_summary():
    """
    Generate snmCT pipeline MappingSummary.csv.gz and save into cwd

    Returns
    -------
    pd.DataFrame
    """
    all_stats = []

    # fastq trimming stats
    df = parse_single_stats_set(f'fastq/*.trimmed.stats.txt',
                                cell_parser_cutadapt_trim_stats)
    all_stats.append(df)

    # hisat-3n DNA mapping
    df = parse_single_stats_set(f'bam/*.hisat3n_dna_summary.txt',
                                cell_parser_hisat_summary, prefix='DNA')
    all_stats.append(df)

    # hisat-3n RNA mapping
    df = parse_single_stats_set(f'rna_bam/*.hisat3n_rna_summary.txt',
                                cell_parser_hisat_summary, prefix='RNA')
    all_stats.append(df)

    # uniquely mapped reads dedup
    df = parse_single_stats_set(f'bam/*.unique_align.deduped.matrix.txt',
                                cell_parser_picard_dedup_stat, prefix='DNAUniqueAlign')
    all_stats.append(df)

    # multi mapped reads dedup
    df = parse_single_stats_set(f'bam/*.multi_align.deduped.matrix.txt',
                                cell_parser_picard_dedup_stat, prefix='DNAMultiAlign')
    all_stats.append(df)

    # uniquely mapped dna reads selection
    df = parse_single_stats_set('bam/*.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv',
                                cell_parser_reads_mc_frac_profile, prefix='UniqueAlign')
    all_stats.append(df)

    # multi mapped dna reads selection
    df = parse_single_stats_set('bam/*.hisat3n_dna.multi_align.deduped.dna_reads.reads_mch_frac.csv',
                                cell_parser_reads_mc_frac_profile, prefix='MultiAlign')
    all_stats.append(df)

    # uniquely mapped rna reads selection
    df = parse_single_stats_set('rna_bam/*.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv',
                                cell_parser_reads_mc_frac_profile)
    all_stats.append(df)

    # allc count
    df = parse_single_stats_set(f'allc/*.allc.tsv.gz.count.csv',
                                cell_parser_allc_count)
    all_stats.append(df)

    # feature count
    df = parse_single_stats_set(f'rna_bam/*.feature_count.tsv.summary',
                                cell_parser_feature_count_summary)
    all_stats.append(df)

    # concatenate all stats
    all_stats = pd.concat(all_stats, axis=1)
    all_stats.index.name = 'cell'
    all_stats.to_csv(f'MappingSummary.csv.gz')
    return all_stats


def snm3c_summary():
    """
    Generate snm3C pipeline MappingSummary.csv.gz and save into cwd

    Returns
    -------
    pd.DataFrame
    """
    all_stats = []

    # fastq trimming stats
    df = parse_single_stats_set(f'fastq/*.trimmed.stats.txt',
                                cell_parser_cutadapt_trim_stats)
    all_stats.append(df)

    # hisat-3n mapping PE
    df = parse_single_stats_set(f'bam/*.hisat3n_dna_summary.txt',
                                cell_parser_hisat_summary)
    all_stats.append(df)

    # hisat-3n mapping split-reads SE
    df = parse_single_stats_set(f'bam/*.hisat3n_dna_split_reads_summary.txt',
                                cell_parser_hisat_summary, prefix='SplitReads')
    all_stats.append(df)

    # uniquely mapped reads dedup
    df = parse_single_stats_set(f'bam/*.all_reads.deduped.matrix.txt',
                                cell_parser_picard_dedup_stat, prefix='UniqueAlign')
    all_stats.append(df)

    # call chromatin contacts
    df = parse_single_stats_set(f'hic/*.all_reads.contact_stats.csv',
                                cell_parser_call_chromatin_contacts)
    all_stats.append(df)

    # allc count
    df = parse_single_stats_set(f'allc/*.allc.tsv.gz.count.csv',
                                cell_parser_allc_count)
    all_stats.append(df)

    # concatenate all stats
    all_stats = pd.concat(all_stats, axis=1)
    all_stats.index.name = 'cell'
    all_stats.to_csv(f'MappingSummary.csv.gz')
    return all_stats
