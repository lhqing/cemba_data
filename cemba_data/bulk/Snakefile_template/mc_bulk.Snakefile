
# Example (required) parameters
# merge_allc_cpu = 10
# mcg_context = 'CGN'
# mch_context = 'CHN'
# bigwig_mch_bin_size = 50
# bigwig_mcg_bin_size = 1
# chrom_size_path = 'PATH_TO_CHROM_SIZE_FILE'
# group = 'GROUP_NAME'

# the main rule is the final target
rule main:
    input:
        f"{group}.{mcg_context}-both.frac.bw",
        f"{group}.{mcg_context}-both.cov.bw",
        f"{group}.{mch_context}-both.frac.bw",
        f"{group}.{mch_context}-both.cov.bw",
        f"{group}.{mcg_context}-Merge.allc.tsv.gz"


# Merge ALLC
rule merge_allc:
    input:
        f"{group}.allc_paths.txt"
    output:
        allc=f"{group}.allc.tsv.gz",
        tbi=f"{group}.allc.tsv.gz.tbi"
    threads:
        max(1, min(int(1.1 * merge_allc_cpu), int(workflow.cores / 1.1)))
    resources:
        mem_mb=merge_allc_cpu * 5000
    shell:
        "allcools merge-allc "
        "--allc_paths {input} "
        "--output_path {output.allc} "
        "--chrom_size_path {chrom_size_path} "
        "--cpu {threads}"


# Extract mCG ALLC for DMR calling
rule extract_allc_mcg:
    input:
        f"{group}.allc.tsv.gz"
    output:
        allc_cg=f"{group}.{mcg_context}-Merge.allc.tsv.gz",
        allc_cg_tbi=f"{group}.{mcg_context}-Merge.allc.tsv.gz.tbi"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools extract-allc "
        "--allc_path {input} "
        "--output_prefix {group} "
        "--mc_contexts {mcg_context} "
        "--chrom_size_path {chrom_size_path} "
        "--strandness merge "
        "--output_format allc "
        "--cpu {threads}"


# Generate mCH BigWig files
rule bigwig_ch:
    input:
        f"{group}.allc.tsv.gz"
    output:
        f"{group}.{mch_context}-both.cov.bw",
        f"{group}.{mch_context}-both.frac.bw"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools allc-to-bigwig "
        "--allc_path {input} "
        "--output_prefix {group} "
        "--bin_size {bigwig_mch_bin_size} "
        "--mc_contexts {mch_context} "
        "--chrom_size_path {chrom_size_path}"


# Generate mCG BigWig files
rule bigwig_cg:
    input:
        f"{group}.allc.tsv.gz"
    output:
        f"{group}.{mcg_context}-both.cov.bw",
        f"{group}.{mcg_context}-both.frac.bw"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools allc-to-bigwig "
        "--allc_path {input} "
        "--output_prefix {group} "
        "--bin_size {bigwig_mcg_bin_size} "
        "--mc_contexts {mcg_context} "
        "--chrom_size_path {chrom_size_path}"
