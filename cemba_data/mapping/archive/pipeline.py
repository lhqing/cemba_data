import pathlib

import pandas as pd
import subprocess
from cemba_data.mapping.archive.bismark_mapping import bismark_mapping
from cemba_data.mapping.archive.generate_allc import generate_allc
from cemba_data.mapping.mct.mct_bismark_bam_filter import dna_reads_selection
from cemba_data.mapping.archive.star_mapping import star_mapping
from cemba_data.utilities import get_configuration


def write_batch_snakefiles(output_dir, mode, cpu=20):
    output_dir = pathlib.Path(output_dir).absolute()
    snakemake_dir = output_dir / 'snakemake'
    n_batches = pd.read_csv(snakemake_dir / 'bismark_bam_list.txt', index_col=0, squeeze=True).max() + 1

    snakemake_commands = []
    for batch_id in range(n_batches):
        if mode == 'mc':
            summary_rule = f"""
include: "{snakemake_dir}/snakefile_generate_allc_{batch_id}"
rule batch_{batch_id}:
    input:
        ["{snakemake_dir}/bismark_mapping_done_{batch_id}",
         "{snakemake_dir}/generate_allc_done_{batch_id}"]

"""
        elif mode == 'mct':
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
        else:
            raise

        with open(snakemake_dir / f'snakefile_batch_{batch_id}', 'w') as f:
            f.write(summary_rule)

        command = f"snakemake -d {snakemake_dir}/{batch_id} --snakefile {snakemake_dir}/snakefile_batch_{batch_id} --cores {cpu} -F"
        snakemake_commands.append(command)

    with open(snakemake_dir / 'snakemake_commands.txt', 'w') as f:
        f.write('\n'.join(snakemake_commands))
    return


def pipeline(output_dir, config_path, batch_size=96, cpu=20, unmapped_fastq=False):
    config = get_configuration(config_path)
    mode = config['mode'].lower()

    # basic mapping
    bismark_mapping(
        output_dir,
        bismark_reference=config['bismark_reference'],
        read_min=int(config['read_min']),
        read_max=int(config['read_max']),
        r1_adapter=config['r1_adapter'],
        r1_left_cut=int(config['r1_left_cut']),
        r1_right_cut=int(config['r1_right_cut']),
        r2_adapter=config['r2_adapter'],
        r2_left_cut=int(config['r2_left_cut']),
        r2_right_cut=int(config['r2_right_cut']),
        split_batch=True,
        batch_size=batch_size,
        unmapped_fastq=unmapped_fastq,
        demultiplex_result=None,  # this is for internal use
        batch_id=None,  # this is for internal use
        testing=None)

    if config['mode'] == 'mct':
        dna_reads_selection(output_dir,
                            mc_rate_max_threshold=float(
                                config['mc_rate_max_threshold']),
                            cov_min_threshold=int(config['cov_min_threshold']))

    generate_allc(output_dir=output_dir,
                  reference_fasta=config['reference_fasta'],
                  num_upstr_bases=int(config['num_upstr_bases']),
                  num_downstr_bases=int(config['num_downstr_bases']),
                  compress_level=int(config['compress_level']))

    if config['mode'] == 'mct':
        star_mapping(output_dir,
                     star_reference=config['star_reference'],
                     gtf_path=config['gtf_path'],
                     threads=cpu,
                     mc_rate_min_threshold=float(config['mc_rate_min_threshold']),
                     cov_min_threshold=int(config['cov_min_threshold']),
                     feature_type=config['feature_type'],
                     id_type=config['id_type'])

    # make a copy of mapping config
    subprocess.run(['cp', str(config_path), f'{output_dir}/snakemake/mapping_config.ini'], check=True)

    write_batch_snakefiles(output_dir=output_dir, mode=mode, cpu=cpu)
    return
