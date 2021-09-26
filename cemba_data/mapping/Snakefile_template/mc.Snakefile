
# Snakemake rules below
# suitable for snmC-seq2, snmC-seq3, NOMe-seq

# use diff mcg_context for normal mC or NOMe
mcg_context = 'CGN' if num_upstr_bases == 0 else 'HCGN'

# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        expand("allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz", cell_id=CELL_IDS,
               mcg_context=mcg_context),
        # also add all the stats path here,
        # once summary is generated, snakemake will delete these stats
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R1.trimmed.stats.tsv", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R2.trimmed.stats.tsv", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS),
    output:
        "MappingSummary.csv.gz"
    shell:
        "yap-internal summary --output_dir ./"

# Trim reads
rule trim_r1:
    input:
        "fastq/{cell_id}-R1.fq.gz"
    output:
        fq=temp("fastq/{cell_id}-R1.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}-R1.trimmed.stats.tsv")
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r1_adapter} {input} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r1_left_cut} -u -{r1_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

rule trim_r2:
    input:
        "fastq/{cell_id}-R2.fq.gz"
    output:
        fq=temp("fastq/{cell_id}-R2.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}-R2.trimmed.stats.tsv")
    threads:
        2
    shell:
        "cutadapt --report=minimal -a {r2_adapter} {input} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r2_left_cut} -u -{r2_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

# bismark mapping, R1 and R2 separately
rule bismark_r1:
    input:
        "fastq/{cell_id}-R1.trimmed.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R1.trimmed_bismark_bt2.bam"),
        stats=temp("bam/{cell_id}-R1.trimmed_bismark_bt2_SE_report.txt")
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R1 with --pbat mode
        "bismark {bismark_reference} {unmapped_param_str} --bowtie2 {input} "
        "--pbat -o bam/ --temp_dir bam/"

rule bismark_r2:
    input:
        "fastq/{cell_id}-R2.trimmed.fq.gz"
    output:
        bam=temp("bam/{cell_id}-R2.trimmed_bismark_bt2.bam"),
        stats=temp("bam/{cell_id}-R2.trimmed_bismark_bt2_SE_report.txt")
    threads:
        3
    resources:
        mem_mb=14000
    shell:
        # map R2 with normal SE mode
        "bismark {bismark_reference} {unmapped_param_str} --bowtie2 {input} "
        "-o bam/ --temp_dir bam/"

# filter bam
rule filter_r1_bam:
    input:
        "bam/{cell_id}-R1.trimmed_bismark_bt2.bam"
    output:
        temp("bam/{cell_id}-R1.trimmed_bismark_bt2.filter.bam")
    shell:
        "samtools view -b -h -q 10 -o {output} {input}"

rule filter_r2_bam:
    input:
        "bam/{cell_id}-R2.trimmed_bismark_bt2.bam"
    output:
        temp("bam/{cell_id}-R2.trimmed_bismark_bt2.filter.bam")
    shell:
        "samtools view -b -h -q 10 -o {output} {input}"

# sort bam
rule sort_r1_bam:
    input:
        "bam/{cell_id}-R1.trimmed_bismark_bt2.filter.bam"
    output:
        temp("bam/{cell_id}-R1.trimmed_bismark_bt2.sorted.bam")
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input}"

rule sort_r2_bam:
    input:
        "bam/{cell_id}-R2.trimmed_bismark_bt2.filter.bam"
    output:
        temp("bam/{cell_id}-R2.trimmed_bismark_bt2.sorted.bam")
    resources:
        mem_mb=1000
    shell:
        "samtools sort -o {output} {input}"

# remove PCR duplicates
rule dedup_r1_bam:
    input:
        "bam/{cell_id}-R1.trimmed_bismark_bt2.sorted.bam"
    output:
        bam=temp("bam/{cell_id}-R1.trimmed_bismark_bt2.deduped.bam"),
        stats=temp("bam/{cell_id}-R1.trimmed_bismark_bt2.deduped.matrix.txt")
    resources:
        mem_mb=1000
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"

rule dedup_r2_bam:
    input:
        "bam/{cell_id}-R2.trimmed_bismark_bt2.sorted.bam"
    output:
        bam=temp("bam/{cell_id}-R2.trimmed_bismark_bt2.deduped.bam"),
        stats=temp("bam/{cell_id}-R2.trimmed_bismark_bt2.deduped.matrix.txt")
    resources:
        mem_mb=1000
    shell:
        "picard MarkDuplicates I={input} O={output.bam} M={output.stats} "
        "REMOVE_DUPLICATES=true TMP_DIR=bam/temp/"

# merge R1 and R2, get final bam
rule merge_bam:
    input:
        "bam/{cell_id}-R1.trimmed_bismark_bt2.deduped.bam",
        "bam/{cell_id}-R2.trimmed_bismark_bt2.deduped.bam"
    output:
        "bam/{cell_id}.final.bam"
    shell:
        "samtools merge -f {output} {input}"

# generate ALLC
rule allc:
    input:
        "bam/{cell_id}.final.bam"
    output:
        allc="allc/{cell_id}.allc.tsv.gz",
        stats=temp("allc/{cell_id}.allc.tsv.gz.count.csv")
    threads:
        2
    resources:
        mem_mb=500
    shell:
        'allcools bam-to-allc '
        '--bam_path {input} '
        '--reference_fasta {reference_fasta} '
        '--output_path {output.allc} '
        '--cpu 1 '
        '--num_upstr_bases {num_upstr_bases} '
        '--num_downstr_bases {num_downstr_bases} '
        '--compress_level {compress_level} '
        '--save_count_df'


# CGN extraction from ALLC
rule cgn_extraction:
    input:
        "allc/{cell_id}.allc.tsv.gz",
    output:
        "allc-{mcg_context}/{cell_id}.{mcg_context}-Merge.allc.tsv.gz",
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
        '--chrom_size_path {chrom_sizes_file} '
