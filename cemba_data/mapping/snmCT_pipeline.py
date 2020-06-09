from .demultiplex import demultiplex_pipeline
from .bismark_mapping import bismark_mapping
from .generate_allc import generate_allc
from .mct_bismark_bam_filter import dna_reads_selection
from .star_mapping import star_mapping, summary_rna_mapping
from .mct_bismark_bam_filter import summarize_select_dna_reads
from .utilities import snakemake, get_configuration

import pathlib
import pandas as pd


def write_batch_snakefiles_ecker_lab(output_dir, cpu=36):
    output_dir = pathlib.Path(output_dir).absolute()
    snakemake_dir = output_dir / 'snakemake'
    n_batches = pd.read_csv(snakemake_dir / 'bismark_bam_list.txt', index_col=0, squeeze=True).max() + 1
    for batch_id in range(n_batches):
        summary_rule = f"""
include: "{snakemake_dir}/snakefile_generate_allc_{batch_id}"
include: "{snakemake_dir}/snakefile_star_mapping_{batch_id}"

rule batch_{batch_id}:
    input:
        ["{snakemake_dir}/bismark_mapping_done_{batch_id}",
         "{snakemake_dir}/select_dna_done_{batch_id}",
         "{snakemake_dir}/generate_allc_done_{batch_id}",
         "{snakemake_dir}/star_mapping_done_{batch_id}"]

"""
        with open(snakemake_dir / f'snakefile_batch_{batch_id}', 'w') as f:
            f.write(summary_rule)

        qsub_temp = f"""
#!/bin/bash
#$ -N batch_{batch_id}
#$ -V 
#$ -l h_rt=999:99:99
#$ -l s_rt=999:99:99
#$ -wd {snakemake_dir}
#$ -e {snakemake_dir}/batch_{batch_id}.error.log
#$ -o {snakemake_dir}/batch_{batch_id}.output.log
#$ -pe smp {cpu}
#$ -l h_vmem=5G

echo JOB_CMD_START batch_{batch_id} COMMAND 0 $(date +"%H:%M:%S-%D")
snakemake -d {snakemake_dir}/{batch_id} --snakefile {snakemake_dir}/snakefile_batch_{batch_id} --cores {cpu}
echo JOB_CMD_RETURN_CODE batch_{batch_id} COMMAND 0 $?
echo JOB_CMD_END batch_{batch_id} COMMAND 0 $(date +"%H:%M:%S-%D")
"""
        with open(snakemake_dir / f'qsub_batch_{batch_id}.sh', 'w') as f:
            f.write(qsub_temp)
    return
