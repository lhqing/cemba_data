
# Example (required) parameters
# merge_allc_cpu = 10
# extract_mcg = False
# mcg_context = 'CGN'
# bigwig_contexts = ['CHN', 'CGN']
# bigwig_bin_size = 50
# chrom_size_path = 'PATH_TO_CHROM_SIZE_FILE'
# groups = [GROUP NAME LIST FOR ALLC MERGE]

# the main rule is the final target
# I only put bigwig files here, so the snakemake can easily handle partial run by using
# --batch main=1/n
rule main:
    input:
        expand("bw/{group}.{mc_context}-Both.{track_type}.bw", group=groups,
               mc_context=bigwig_contexts, track_type=['cov', 'rate'])

# Merge ALLC
rule merge_allc:
    input:
        "allc/{group}.allc_paths.txt"
    output:
        allc="allc/{group}.allc.tsv.gz",
        tbi="allc/{group}.allc.tsv.gz.tbi"
    threads:
        max(1, min(int(1.5 * merge_allc_cpu), int(workflow.cores / 1.5)))
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
if extract_mcg:
    rule extract_allc_mcg:
        input:
            "allc/{group}.allc.tsv.gz"
        output:
            allc_cg="allc/{group}.{mcg_context}-Merge.allc.tsv.gz",
            allc_cg_tbi="allc/{group}.{mcg_context}-Merge.allc.tsv.gz.tbi"
        params:
            prefix="allc/{group}"
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
    def bigwig_input(wildcards):
        input = {'allc': f"allc/{wildcards.group}.allc.tsv.gz",
                 'allc_cg': f"allc/{wildcards.group}.{mcg_context}-Merge.allc.tsv.gz"}
        return input
else:
    def bigwig_input(wildcards):
        input = {'allc': f"allc/{wildcards.group}.allc.tsv.gz"}
        return input


# Generate BigWig for browsers
# Note that we also put the extract_mcg output as input here, not because they are needed
# but just for the rule dependency so that the main rule is aware of the extract_mcg rules
rule bigwig:
    input:
        unpack(bigwig_input)
    output:
        ["bw/{group}.{mc_context}-Both.cov.bw"
         for mc_context in bigwig_contexts],
        ["bw/{group}.{mc_context}-Both.rate.bw"
         for mc_context in bigwig_contexts],
    params:
        prefix="bw/{group}",
        context_str=' '.join(bigwig_contexts)
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools allc-to-bigwig "
        "--allc_path {input.allc} "
        "--bin_size {bigwig_bin_size} "
        "--output_prefix {params.prefix} "
        "--chrom_size_path {chrom_size_path} "
        "--mc_contexts {params.context_str} "
        "--remove_additional_chrom "
        "--cpu {threads}"
