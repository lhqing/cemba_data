
# Snakemake rules below
# suitable for snmCT-seq, snmC2T-seq (with NOMe treatment)

# the summary rule is the final target
rule summary:
    input:
        expand("allc/{cell_id}.allc.tsv.gz", cell_id=CELL_IDS),
        # also add all the stats path here,
        # once summary is generated, snakemake will delete these stats
        expand("allc/{cell_id}.allc.tsv.gz.count.csv", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R1.trimmed.stats.txt", cell_id=CELL_IDS),
        expand("fastq/{cell_id}-R2.trimmed.stats.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.trimmed_bismark_bt2.deduped.matrix.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R1.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}-R2.trimmed_bismark_bt2_SE_report.txt", cell_id=CELL_IDS),
        expand("bam/{cell_id}.dna_reads.bam.reads_profile.csv", cell_id=CELL_IDS),
        'rna_bam/TotalRNAAligned.filtered.bam',  # needed for count star mapped reads by RG
        'rna_bam/TotalRNALog.final.out',
        'rna_bam/TotalRNALog.out',
        'rna_bam/TotalRNALog.progress.out',
        'rna_bam/TotalRNAAligned.rna_reads.bam',
        'rna_bam/TotalRNAAligned.rna_reads.bam.reads_profile.csv',
        'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv',
        'rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary'
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
        stats=temp("fastq/{cell_id}-R1.trimmed.stats.txt")
    threads:
        2
    shell:
        "cutadapt -a R1Adapter={r1_adapter} "
        "-a TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG "
        "-a N6=AAGCAGTGGTATCAACGCAGAGTAC "
        "-a TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT "
        "-a N6_rc=GTACTCTGCGTTGATACCACTGCTT "
        "-a 3PpolyT=TTTTTTTTTTTTTTTX "
        "-g 5PpolyT=XTTTTTTTTTTTTTTT "
        "-a 3PpolyA=AAAAAAAAAAAAAAAX "
        "-g 5PpolyA=XAAAAAAAAAAAAAAA "
        "-a polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "
        "-a polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
        "-a ISPCR_F=AAGCAGTGGTATCAACGCAGAGT "
        "-a ISPCR_R=ACTCTGCGTTGATACCACTGCTT "
        "{input} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r1_left_cut} -u -{r1_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

rule trim_r2:
    input:
        "fastq/{cell_id}-R2.fq.gz"
    output:
        fq=temp("fastq/{cell_id}-R2.trimmed.fq.gz"),
        stats=temp("fastq/{cell_id}-R2.trimmed.stats.txt")
    threads:
        2
    shell:
        "cutadapt -a R2Adapter={r2_adapter} "
        "-a TSO=AAGCAGTGGTATCAACGCAGAGTGAATGG "
        "-a N6=AAGCAGTGGTATCAACGCAGAGTAC "
        "-a TSO_rc=CCATTCACTCTGCGTTGATACCACTGCTT "
        "-a N6_rc=GTACTCTGCGTTGATACCACTGCTT "
        "-a 3PpolyT=TTTTTTTTTTTTTTTX "
        "-g 5PpolyT=XTTTTTTTTTTTTTTT "
        "-a 3PpolyA=AAAAAAAAAAAAAAAX "
        "-g 5PpolyA=XAAAAAAAAAAAAAAA "
        "-a polyTLong=TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT "
        "-a polyALong=AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA "
        "-a ISPCR_F=AAGCAGTGGTATCAACGCAGAGT "
        "-a ISPCR_R=ACTCTGCGTTGATACCACTGCTT "
        "{input} 2> {output.stats} | "
        "cutadapt --report=minimal -O 6 -q 20 -u {r2_left_cut} -u -{r2_right_cut} -m 30 "
        "-o {output.fq} - >> {output.stats}"

# first is bismark mapping, R1 and R2 separately
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
        temp("bam/{cell_id}.final.bam")
    shell:
        "samtools merge -f {output} {input}"

# select DNA reads from final bismark mapped bam
rule select_dna:
    input:
        "bam/{cell_id}.final.bam"
    output:
        bam="bam/{cell_id}.dna_reads.bam",
        stats=temp('bam/{cell_id}.dna_reads.bam.reads_profile.csv')
    shell:
        'yap-internal select-dna-reads --input_bam {input} '
        '--output_bam {output.bam} --mc_rate_max_threshold {mc_rate_max_threshold} '
        '--cov_min_threshold {dna_cov_min_threshold} '
        '{nome_flag_str} '
        '--assay_type mc'

# generate ALLC using dna_reads.bam
rule allc:
    input:
        "bam/{cell_id}.dna_reads.bam"
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


# RNA mapping, also start from trimmed fastq
cell_ids_str = ' , ID:'.join(CELL_IDS)
# star separate multiple input by ,
star_input_str = ','.join([f"fastq/{cell_id}-R1.trimmed.fq.gz" for cell_id in CELL_IDS])

rule star:
    input:
        # here we only use R1 SE for RNA,
        # R2 SE or R1R2 PE is worse than R1 actually, due to R2's low quality
        # And we map all cells together, so the genome is only load once
        # each cell will have a different @RG tag
        expand("fastq/{cell_id}-R1.trimmed.fq.gz", cell_id = CELL_IDS)
    output:
        temp('rna_bam/TotalRNAAligned.out.bam'),
        temp('rna_bam/TotalRNALog.final.out'),
        temp('rna_bam/TotalRNALog.out'),
        temp('rna_bam/TotalRNALog.progress.out'),
        temp('rna_bam/TotalRNASJ.out.tab')
    threads:
        workflow.cores * 0.8  # workflow.cores is user provided cores for snakemake
    resources:
        mem_mb=48000
    shell:
        'STAR --runThreadN {threads} '
        '--genomeDir {star_reference} '
        '--alignEndsType Local '
        '--genomeLoad NoSharedMemory '
        '--outSAMstrandField intronMotif '
        '--outSAMtype BAM Unsorted '
        '--outSAMunmapped None '
        '--outSAMattributes NH HI AS NM MD '
        '--sjdbOverhang 100 '
        '--outFilterType BySJout '  # ENCODE standard options
        '--outFilterMultimapNmax 20 '  # ENCODE standard options
        '--alignSJoverhangMin 8 '  # ENCODE standard options
        '--alignSJDBoverhangMin 1 '  # ENCODE standard options
        '--outFilterMismatchNmax 999 '  # ENCODE standard options
        '--outFilterMismatchNoverLmax 0.04 '  # ENCODE standard options
        '--alignIntronMin 20 '  # ENCODE standard options
        '--alignIntronMax 1000000 '  # ENCODE standard options
        '--alignMatesGapMax 1000000 '  # ENCODE standard options
        '--outFileNamePrefix rna_bam/TotalRNA '
        '--readFilesIn {star_input_str} '
        '--readFilesCommand gzip -cd '
        '--outSAMattrRGline ID:{cell_ids_str}'

rule filter_bam:
    input:
        'rna_bam/TotalRNAAligned.out.bam'
    output:
        temp('rna_bam/TotalRNAAligned.filtered.bam')
    threads:
        min(workflow.cores * 0.8, 10)
    shell:
        "samtools sort -@ {threads} -m 2G {input} | samtools view -bh -q 10 -o {output} -"

rule select_rna:
    input:
        'rna_bam/TotalRNAAligned.filtered.bam'
    output:
        bam='rna_bam/TotalRNAAligned.rna_reads.bam',
        stats=temp('rna_bam/TotalRNAAligned.rna_reads.bam.reads_profile.csv')
    shell:
        'yap-internal select-rna-reads ' \
        '--input_bam {input} ' \
        '--output_bam {output.bam} ' \
        '--mc_rate_min_threshold {mc_rate_min_threshold} ' \
        '--cov_min_threshold {rna_cov_min_threshold} '
        '{nome_flag_str} '
        '--assay_type mc'

rule feature_count:
    input:
        'rna_bam/TotalRNAAligned.rna_reads.bam'
    output:
        count='rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv',
        stats=temp('rna_bam/TotalRNAAligned.rna_reads.feature_count.tsv.summary')
    threads:
        min(workflow.cores * 0.8, 10)
    resources:
        mem_mb=1000
    shell:
        'featureCounts -t {feature_type} -g {id_type} ' \
        '-a {gtf_path} -o {output.count} --byReadGroup -T {threads} {input}'
