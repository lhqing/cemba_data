# ==================================================
# Import
# ==================================================


import yaml
import pathlib
from cemba_data.hisat3n import *


# ==================================================
# Preparation
# ==================================================


# read mapping config and put all variables into the locals()
config, config_dict = read_mapping_config()
# print('Usings these mapping parameters:')
# for _k, _v in config_dict.items():
#     print(f'{_k} = {_v}')

# fastq table and cell IDs
fastq_table = validate_cwd_fastq_paths()
CELL_IDS = fastq_table.index.tolist()

# print(f"Found {len(CELL_IDS)} FASTQ pairs in fastq/ .")

mcg_context = 'CGN' if int(config.num_upstr_bases) == 0 else 'HCGN'
repeat_index_flag = "--repeat" if config.hisat3n_repeat_index_type == 'repeat' else "--no-repeat-index"


# ==================================================
# Mapping summary
# ==================================================


# the summary rule is the final target
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt", cell_id=CELL_IDS),
        # dna mapping
        expand("bam/{cell_id}.hisat3n_dna_summary.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv", cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam.bai", cell_id=CELL_IDS),
        # rna mapping
        expand("rna_bam/{cell_id}.hisat3n_rna_summary.txt", cell_id=CELL_IDS),
        expand("rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv", cell_id=CELL_IDS),
        expand("rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam.bai", cell_id=CELL_IDS),
        expand("rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv", cell_id=CELL_IDS),
        expand("rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv.summary", cell_id=CELL_IDS),
        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",
               cell_id=CELL_IDS, mcg_context=mcg_context),
    output:
        "MappingSummary.csv.gz"
    run:
        snmct_summary()
        aggregate_feature_counts()

        # cleanup
        shell("rm -rf bam/temp")


# ==================================================
# FASTQ Trimming
# ==================================================


# Trim reads
# sort the fastq files so that R1 and R2 are in the same order
rule sort_R1:
    input:
        "fastq/{cell_id}-R1.fq.gz",
    output:
        temp("fastq/{cell_id}-R1_sort.fq")
    threads:
        1.5
    shell:
        'zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output} '


rule sort_R2:
    input:
        "fastq/{cell_id}-R2.fq.gz",
    output:
        temp("fastq/{cell_id}-R2_sort.fq")
    threads:
        1.5
    shell:
        'zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > {output} '


rule trim:
    input:
        R1="fastq/{cell_id}-R1_sort.fq",
        R2="fastq/{cell_id}-R2_sort.fq"
    output:
        R1=temp("fastq/{cell_id}-R1.trimmed.fq.gz"),
        R2=temp("fastq/{cell_id}-R2.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}.trimmed.stats.txt")
    threads:
        1
    shell:
        "cutadapt "
        "-a R1Adapter={config.r1_adapter} "
        "-a TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG "
        "-a TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT "
        "-a N6=AAGCAGTGGTATCAACGCAGAGTAC "
        "-a N6_rc=GTACTCTGCGTTGATACCACTGCTT "
        "-a 3PpolyT=TTTTTTTTTTTTTTTX "
        "-a 3PpolyA=AAAAAAAAAAAAAAAX "
        "-a polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "
        "-a polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
        "-a ISPCR_F=AAGCAGTGGTATCAACGCAGAGT "
        "-a ISPCR_R=ACTCTGCGTTGATACCACTGCTT "
        "-A R2Adapter={config.r2_adapter} "
        "-A TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG "
        "-A TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT "
        "-A N6=AAGCAGTGGTATCAACGCAGAGTAC "
        "-A N6_rc=GTACTCTGCGTTGATACCACTGCTT "
        "-A 3PpolyT=TTTTTTTTTTTTTTTX "
        "-A 3PpolyA=AAAAAAAAAAAAAAAX "
        "-A polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "
        "-A polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
        "-A ISPCR_F=AAGCAGTGGTATCAACGCAGAGT "
        "-A ISPCR_R=ACTCTGCGTTGATACCACTGCTT "
        "-g 5PpolyT=XTTTTTTTTTTTTTTT "
        "-g 5PpolyA=XAAAAAAAAAAAAAAA "
        "-G 5PpolyT=XTTTTTTTTTTTTTTT "
        "-G 5PpolyA=XAAAAAAAAAAAAAAA "
        "--report=minimal "
        "-O 6 "
        "-q 20 "
        "-u {config.r1_left_cut} "
        "-u -{config.r1_right_cut} "
        "-U {config.r2_left_cut} "
        "-U -{config.r2_right_cut} "
        "-m 30:30 "
        "--pair-filter 'both' "
        "-o {output.R1} "
        "-p {output.R2} "
        "{input.R1} {input.R2} "
        "> {output.stats}"


# ==================================================
# HISAT-3N DNA Mapping
# ==================================================


# Paired-end Hisat3n mapping using DNA mode
rule hisat_3n_pairend_mapping_dna_mode:
    input:
        R1="fastq/{cell_id}-R1.trimmed.fq.gz",
        R2="fastq/{cell_id}-R2.trimmed.fq.gz"
    output:
        bam=temp("bam/{cell_id}.hisat3n_dna.unsort.bam"),
        stats=temp("bam/{cell_id}.hisat3n_dna_summary.txt")
    threads:
        8
    resources:
        mem_mb=8000
    shell:
        "hisat-3n "
        "{config.hisat3n_dna_reference} "
        "-q "
        "-1 {input.R1} "
        "-2 {input.R2} "
        # "--directional-mapping "  # this can speed up 2X as the snmC reads are directional
        "--base-change C,T "
        "{repeat_index_flag} "
        "--no-spliced-alignment "  # this is important for DNA mapping
        "--no-temp-splicesite "
        "-t "
        "--new-summary "
        "--summary-file {output.stats} "
        "--threads {threads} "
        "| "
        "samtools view "
        "-b -q 10 -o {output.bam}"  # -q 10 will filter out multi-aligned reads


rule sort_dna_bam:
    input:
        "bam/{cell_id}.hisat3n_dna.unsort.bam"
    output:
        temp("bam/{cell_id}.hisat3n_dna.unique_align.bam")
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        "samtools sort -O BAM -o {output} {input}"


# remove PCR duplicates
rule dedup_unique_bam:
    input:
        "bam/{cell_id}.hisat3n_dna.unique_align.bam"
    output:
        bam=temp("bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam"),
        stats=temp("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt")
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"


rule select_unique_bam_dna_reads:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam"
    output:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam",
        stats=temp("bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.reads_mch_frac.csv")
    resources:
        mem_mb=100
    run:
        select_mct_reads_normal(
            input_bam=input.bam,
            output_bam=output.bam,
            mode='dna',
            mc_rate_max_threshold=0.5,
            cov_min_threshold=3,
            nome=False
        )


rule index_unique_bam_dna_reads:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam"
    output:
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam.bai"
    shell:
        "samtools index {input.bam}"


# ==================================================
# HISAT-3N RNA Mapping
# ==================================================


# Paired-end Hisat3n mapping using RNA mode
rule hisat_3n_pairend_mapping_rna_mode:
    input:
        R1="fastq/{cell_id}-R1.trimmed.fq.gz",
        R2="fastq/{cell_id}-R2.trimmed.fq.gz"
    output:
        bam="rna_bam/{cell_id}.hisat3n_rna.unsort.bam",
        stats=temp("rna_bam/{cell_id}.hisat3n_rna_summary.txt")
    threads:
        8
    resources:
        mem_mb=8000
    shell:
        "hisat-3n "
        "-x {config.hisat3n_rna_reference} "
        "-q "
        "-1 {input.R1} "
        "-2 {input.R2} "
        "--base-change C,T "
        "{repeat_index_flag} "
        "-t "
        "--new-summary "
        "--summary-file {output.stats} "
        "--threads {threads} "
        "| "
        "samtools addreplacerg "  # add read group @RG to the reads in order to use featuerCounts
        "-r '@RG\tID:{wildcards.cell_id}' -u -o - -"
        "| "
        "samtools view "
        "-b -q 10 -o {output.bam}"  # -q 10 will filter out multi-aligned reads


rule sort_rna_bam:
    input:
        "rna_bam/{cell_id}.hisat3n_rna.unsort.bam"
    output:
        temp("rna_bam/{cell_id}.hisat3n_rna.bam")
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        "samtools sort -O BAM -o {output} {input}"


# skip dedup step for RNA reads


rule select_unique_bam_rna_reads:
    input:
        bam="rna_bam/{cell_id}.hisat3n_rna.bam"
    output:
        bam="rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam",
        stats=temp("rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.reads_mch_frac.csv")
    resources:
        mem_mb=100
    run:
        select_mct_reads_normal(
            input_bam=input.bam,
            output_bam=output.bam,
            mode='rna',
            mc_rate_min_threshold=0.9,
            cov_min_threshold=3,
            nome=False
        )


rule index_unique_bam_rna_reads:
    input:
        bam="rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam"
    output:
        bai="rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam.bai"
    shell:
        "samtools index {input.bam}"


rule feature_count:
    input:
        'rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.bam'
    output:
        tsv=temp('rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv'),
        stats=temp('rna_bam/{cell_id}.hisat3n_rna.unique_align.rna_reads.feature_count.tsv.summary')
    threads:
        1
    resources:
        mem_mb=1000
    shell:
        'featureCounts -t {config.feature_type} -g {config.id_type} ' \
        '-a {config.gtf_path} -o {output.tsv} --byReadGroup -T {threads} {input}'



# ==================================================
# Generate ALLC
# ==================================================


# generate ALLC
rule unique_reads_allc:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.dna_reads.bam"
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        stats=temp("allc/{cell_id}.allc.tsv.gz.count.csv")
    threads:
        1.5
    resources:
        mem_mb=500
    shell:
        'allcools bam-to-allc '
        '--bam_path {input.bam} '
        '--reference_fasta {config.reference_fasta} '
        '--output_path {output.allc} '
        '--cpu {threads} '
        '--num_upstr_bases {config.num_upstr_bases} '
        '--num_downstr_bases {config.num_downstr_bases} '
        '--compress_level {config.compress_level} '
        '--save_count_df '
        '--convert_bam_strandness '


# CGN extraction from ALLC
rule unique_reads_cgn_extraction:
    input:
        "allc/{cell_id}.allc.tsv.gz",
    output:
        allc="allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",
        tbi="allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",
    params:
        prefix="allc-{mcg_context}/{cell_id}",
    threads:
        1
    resources:
        mem_mb=100
    shell:
        'allcools extract-allc '
        '--strandness merge '
        '--allc_path  {input} '
        '--output_prefix {params.prefix} '
        '--mc_contexts {mcg_context} '
        '--chrom_size_path {config.chrom_size_path} '
