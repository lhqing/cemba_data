import pathlib

from cemba_data.mapping import \
    prepare_select_rna_reads, \
    command_runner, \
    summarize_select_rna_reads, \
    star_mapping, \
    summarize_star_mapping, star_bam_qc


def pipeline_mc_info(output_dir, config_path, qsub=True, cpu=10):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'star_bam'
    bam_dir.mkdir(exist_ok=True, parents=True)

    # prepare STAR mapping
    star_mapping_records, star_mapping_commands = star_mapping(input_dir=fastq_dir,
                                                               output_dir=bam_dir,
                                                               config=config_path)
    # runner STAR mapping
    if qsub:
        raise NotImplementedError
    else:
        command_runner(star_mapping_commands, runner=None, cpu=cpu)

    # summary STAR mapping
    summarize_star_mapping(bam_dir=bam_dir)

    # prepare STAR BAM QC
    star_bam_qc_records, star_bam_qc_commands = star_bam_qc(output_dir=output_dir, config=config_path)

    # runner STAR BAM QC
    if qsub:
        raise NotImplementedError
    else:
        command_runner(star_bam_qc_commands, runner=None, cpu=cpu)

    # prepare select RNA
    select_rna_records, select_rna_commands = prepare_select_rna_reads(output_dir=bam_dir, config=config_path)

    # run select RNA
    if qsub:
        raise NotImplementedError
    else:
        command_runner(select_rna_commands, runner=None, cpu=cpu)

    # summarize select RNA
    summarize_select_rna_reads(output_dir=bam_dir)

    # feature count
    

    # final summary
