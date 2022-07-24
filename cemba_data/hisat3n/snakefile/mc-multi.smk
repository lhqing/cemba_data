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
DEFAULT_CONFIG = {
    'hisat3n_repeat_index_type': '',
    'r1_adapter': 'AGATCGGAAGAGCACACGTCTGAAC',
    'r2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGA',
    'r1_right_cut': 10,
    'r2_right_cut': 10,
    'r1_left_cut': 10,
    'r2_left_cut': 10,
    'min_read_length': 30,
    'num_upstr_bases': 0,
    'num_downstr_bases': 2,
    'compress_level': 5,
    'hisat3n_threads': 11,
    # the post_mapping_script can be used to generate dataset, run other process etc.
    # it gets executed before the final summary function.
    # the default command is just a placeholder that has no effect
    'post_mapping_script': 'true',
}
REQUIRED_CONFIG = ['hisat3n_dna_reference', 'reference_fasta', 'chrom_size_path']

for k, v in DEFAULT_CONFIG.items():
    if k not in config:
        config[k] = v

missing_key = []
for k in REQUIRED_CONFIG:
    if k not in config:
        missing_key.append(k)
if len(missing_key) > 0:
    raise ValueError('Missing required config: {}'.format(missing_key))

# fastq table and cell IDs
fastq_table = validate_cwd_fastq_paths()
CELL_IDS = fastq_table.index.tolist()

mcg_context = 'CGN' if int(config['num_upstr_bases']) == 0 else 'HCGN'
repeat_index_flag = "--repeat" if config['hisat3n_repeat_index_type'] == 'repeat' else "--no-repeat-index"


# ==================================================
# Mapping summary
# ==================================================


# the summary rule is the final target
rule summary:
    input:
        # fastq trim
        expand("fastq/{cell_id}.trimmed.stats.txt",cell_id=CELL_IDS),
        # dna mapping
        expand("bam/{cell_id}.hisat3n_dna_summary.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt",cell_id=CELL_IDS),
        expand("bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt",cell_id=CELL_IDS),
        # allc
        expand("allc/{cell_id}.allc.tsv.gz.count.csv",cell_id=CELL_IDS),
        expand("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv",cell_id=CELL_IDS),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz.tbi",
            cell_id=CELL_IDS,mcg_context=mcg_context),
    output:
        "MappingSummary.csv.gz"
    run:
        # execute any post-mapping script before generating the final summary
        shell(config['post_mapping_script'])

        snmc_summary()

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
        # change to sort_R1 and sort_R2 output if the FASTQ name is disordered
        R1="fastq/{cell_id}-R1.fq.gz",
        R2="fastq/{cell_id}-R2.fq.gz"
    output:
        R1=temp("fastq/{cell_id}-R1.trimmed.fq.gz"),
        R2=temp("fastq/{cell_id}-R2.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}.trimmed.stats.txt")
    threads:
        1
    shell:
        "cutadapt "
        "-a R1Adapter={config[r1_adapter]} "
        "-A R2Adapter={config[r2_adapter]} "
        "--report=minimal "
        "-O 6 "
        "-q 20 "
        "-u {config[r1_left_cut]} "
        "-u -{config[r1_right_cut]} "
        "-U {config[r2_left_cut]} "
        "-U -{config[r2_right_cut]} "
        "-Z "
        "-m {config[min_read_length]}:30 "
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
        "{config[hisat3n_dna_reference]} "
        "-q "
        "-1 {input.R1} "
        "-2 {input.R2} "
        "--directional-mapping-reverse "  # this can speed up 2X as the snmC reads are directional
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
        "-b -q 1 -o {output.bam}"


rule sort_bam:
    input:
        "bam/{cell_id}.hisat3n_dna.unsort.bam"
    output:
        temp("bam/{cell_id}.hisat3n_dna.bam")
    resources:
        mem_mb=1000
    threads:
        1
    shell:
        "samtools sort -O BAM -o {output} {input}"


# Separate unique aligned reads and multi-aligned reads with length > 30
# TODO: make sure how to separate multi-align reads? or reads map to repeat regions in the genome?
# TODO right now, we are just using mapq == 1 as multi-align reads, but this might not be right
rule split_unique_and_multi_align_bam_dna:
    input:
        bam="bam/{cell_id}.hisat3n_dna.bam"
    output:
        unique=temp("bam/{cell_id}.hisat3n_dna.unique_align.bam"),
        multi=temp("bam/{cell_id}.hisat3n_dna.multi_align.bam")
    run:
        separate_unique_and_multi_align_reads(
            in_bam_path=input.bam,
            out_unique_path=output.unique,
            out_multi_path=output.multi,
            out_unmappable_path=None,
            mapq_cutoff=10,
            qlen_cutoff=30
        )


# remove PCR duplicates
rule dedup_unique_bam:
    input:
        "bam/{cell_id}.hisat3n_dna.unique_align.bam"
    output:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        stats=temp("bam/{cell_id}.hisat3n_dna.unique_align.deduped.matrix.txt")
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"


rule dedup_multi_bam:
    input:
        "bam/{cell_id}.hisat3n_dna.multi_align.bam"
    output:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        stats=temp("bam/{cell_id}.hisat3n_dna.multi_align.deduped.matrix.txt")
    resources:
        mem_mb=1000
    threads:
        2
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"


rule index_unique_bam_dna_reads:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam"
    output:
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
    shell:
        "samtools index {input.bam}"


rule index_multi_bam_dna_reads:
    input:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam"
    output:
        bai="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam.bai"
    shell:
        "samtools index {input.bam}"


# ==================================================
# Generate ALLC
# ==================================================


# generate ALLC
rule unique_reads_allc:
    input:
        bam="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.unique_align.deduped.bam.bai"
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
        '--reference_fasta {config[reference_fasta]} '
        '--output_path {output.allc} '
        '--num_upstr_bases {config[num_upstr_bases]} '
        '--num_downstr_bases {config[num_downstr_bases]} '
        '--compress_level {config[compress_level]} '
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
        '--chrom_size_path {config[chrom_size_path]} '


# generate ALLC
rule multi_reads_allc:
    input:
        bam="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam",
        bai="bam/{cell_id}.hisat3n_dna.multi_align.deduped.bam.bai"
    output:
        allc="allc-multi/{cell_id}.allc_multi.tsv.gz",
        stats=temp("allc-multi/{cell_id}.allc_multi.tsv.gz.count.csv")
    threads:
        1.5
    resources:
        mem_mb=500
    shell:
        'allcools bam-to-allc '
        '--bam_path {input.bam} '
        '--reference_fasta {config[reference_fasta]} '
        '--output_path {output.allc} '
        '--num_upstr_bases {config[num_upstr_bases]} '
        '--num_downstr_bases {config[num_downstr_bases]} '
        '--compress_level {config[compress_level]} '
        '--save_count_df '
        '--min_mapq 0 '  # for multi-mapped reads, skip mapq filter
        '--convert_bam_strandness '
