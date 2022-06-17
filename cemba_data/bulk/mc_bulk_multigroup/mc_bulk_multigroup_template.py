MERGE_TEMPLATE = '''
# Example (required) parameters
# merge_allc_cpu = 10
# chrom_size_path = 'PATH_TO_CHROM_SIZE_FILE'
# merge_sample_prefixes = '[]'
# copy_sample_prefixes = '[]'
# group = 'GROUP_NAME'
sample_prefixes = merge_sample_prefixes + copy_sample_prefixes

# the main rule is the final target
rule main:
    input:
        expand("{sample}.allc.tsv.gz", sample=sample_prefixes),
        expand("{sample}.allc.tsv.gz.tbi", sample=sample_prefixes),
#    output:
#        f"{group}.finished"
#    shell:
#        "date > {output}"
        


# Merge ALLC
rule merge_allc:
    input:
        "{sample}.pathlist",
    output:
        allc="{sample}.allc.tsv.gz",
        tbi="{sample}.allc.tsv.gz.tbi"
    threads:
        max(1, min(int(1.1 * merge_allc_cpu), int(workflow.cores / 1.1)))
    resources:
        mem_mb=merge_allc_cpu * 5000
    run:
        if wildcards.sample in merge_sample_prefixes:
            shell("allcools merge-allc "
                  "--allc_paths {input} "
                  "--output_path {output.allc} "
                  "--chrom_size_path {chrom_size_path} "
                  "--cpu {threads}")
        else:
            shell("cp $(cat {input}) {output.allc} ;"
                  "cp $(cat {input}).tbi {output.tbi} ;")

'''

MERGE_EXTRACT_TEMPLATE = '''
# Example (required) parameters
# merge_allc_cpu = 10
# mcg_context = 'CGN'
# chrom_size_path = 'PATH_TO_CHROM_SIZE_FILE'
# merge_sample_prefixes = '[]'
# copy_sample_prefixes = '[]'
# group = 'GROUP_NAME'
sample_prefixes = merge_sample_prefixes + copy_sample_prefixes

# the main rule is the final target
rule main:
    input:
        expand("{sample}.{mcg_context}-Merge.allc.tsv.gz", sample=sample_prefixes, mcg_context=[mcg_context]),
        expand("{sample}.{mcg_context}-Merge.allc.tsv.gz.tbi", sample=sample_prefixes, mcg_context=[mcg_context]),
#    output:
#        f"{group}.finished"
#    shell:
#        "date > {output}"
        

# Merge ALLC
rule merge_allc:
    input:
        "{sample}.pathlist",
    output:
        allc="{sample}.allc.tsv.gz",
        tbi="{sample}.allc.tsv.gz.tbi"
    threads:
        max(1, min(int(1.1 * merge_allc_cpu), int(workflow.cores / 1.1)))
    resources:
        mem_mb=merge_allc_cpu * 5000
    run:
        if wildcards.sample in merge_sample_prefixes:
            shell("allcools merge-allc "
                  "--allc_paths {input} "
                  "--output_path {output.allc} "
                  "--chrom_size_path {chrom_size_path} "
                  "--cpu {threads}")
        else:
            shell("cp $(cat {input}) {output.allc} ;"
                  "cp $(cat {input}).tbi {output.tbi} ;")

# Extract mCG ALLC for DMR calling
rule extract_allc_mcg:
    input:
        "{sample}.allc.tsv.gz"
    output:
        allc_cg="{sample}.{mcg_context}-Merge.allc.tsv.gz",
        allc_cg_tbi="{sample}.{mcg_context}-Merge.allc.tsv.gz.tbi"
    threads:
        1
    resources:
        mem_mb=100
    shell:
        "allcools extract-allc "
        "--allc_path {input} "
        "--output_prefix {wildcards.sample} "
        "--mc_contexts {mcg_context} "
        "--chrom_size_path {chrom_size_path} "
        "--strandness merge "
        "--output_format allc "
        "--cpu {threads}"
'''
