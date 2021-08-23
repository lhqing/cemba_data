
# Example (required) parameters
# merge_allc_cpu = 10
# mcg_context = 'CGN'
# mch_context = 'CHN'
# bigwig_mch_bin_size = 50
# bigwig_mcg_bin_size = 1
# chrom_size_path = 'PATH_TO_CHROM_SIZE_FILE'
# groups = [GROUP NAME LIST FOR ALLC MERGE]

# the main rule is the final target
rule main:
    input:
        expand(f"{{group}}.{mcg_context}-Both.{{track_type}}.bw",
               group=groups, track_type=['cov', 'rate']),
        expand(f"{{group}}.{mch_context}-Both.{{track_type}}.bw",
               group=groups, track_type=['cov', 'rate']),
        expand(f"{{group}}.{mcg_context}-Merge.allc.tsv.gz",
               group=groups)

# Merge ALLC
rule merge_allc:
    input:
        "{group}.allc_paths.txt"
    output:
        allc="{group}.allc.tsv.gz",
        tbi="{group}.allc.tsv.gz.tbi"
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

# Optional: Extract mCG ALLC for DMR calling
# If the ALLC before merging is already extracted, this step is not necessary
rule extract_allc_mcg:
    input:
        "{group}.allc.tsv.gz"
    output:
        allc_cg="{group}.{mcg_context}-Merge.allc.tsv.gz",
        allc_cg_tbi="{group}.{mcg_context}-Merge.allc.tsv.gz.tbi"
    params:
        prefix="{group}"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools extract-allc "
        "--allc_path {input} "
        "--output_prefix {params.prefix} "
        "--mc_contexts {mcg_context} "
        "--chrom_size_path {chrom_size_path} "
        "--strandness merge "
        "--output_format allc "
        "--cpu {threads}"


# Generate mCH BigWig files
rule bigwig_ch:
    input:
        "{group}.allc.tsv.gz"
    output:
        "{group}.{mch_context}-both.cov.bw",
        "{group}.{mch_context}-both.rate.bw",
    params:
        prefix="{group}"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools allc-to-bigwig "
        "--allc_path {input.allc} "
        "--output_prefix {params.prefix} "
        "--bin_size {bigwig_mch_bin_size} "
        "--mc_contexts {mch_context} "
        "--chrom_size_path {chrom_size_path} "
        "--strandness both"

# Generate mCG BigWig files
rule bigwig_cg:
    input:
        "{group}.allc.tsv.gz"
    output:
        "{group}.{mcg_context}-both.cov.bw",
        "{group}.{mcg_context}-both.rate.bw",
    params:
        prefix="{group}"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools allc-to-bigwig "
        "--allc_path {input.allc} "
        "--output_prefix {params.prefix} "
        "--bin_size {bigwig_mcg_bin_size} "
        "--mc_contexts {mcg_context} "
        "--chrom_size_path {chrom_size_path} "
        "--strandness both"
