import pathlib

from cemba_data.mapping import \
    bismark_mapping, \
    command_runner, \
    summarize_bismark_mapping, \
    bismark_bam_qc, \
    summarize_bismark_bam_qc, merge_bam, generate_allc


def pipeline_mc_info(output_dir, config_path, qsub=True):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'bismark_bam'
    bam_dir.mkdir(exist_ok=True, parents=True)

    # prepare bismark mapping
    # bismark_records ['uid', 'index_name', 'read_type', 'bam_path']
    bismark_records, bismark_commands = bismark_mapping(input_dir=fastq_dir,
                                                        output_dir=bam_dir,
                                                        config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        for cmd in bismark_commands:
            command_runner(cmd)

    # summarize bismark mapping
    summarize_bismark_mapping(output_dir=bam_dir)

    # prepare bismark bam qc
    # bam_qc_records ['uid', 'index_name', 'read_type', 'bam_path']
    bam_qc_records, bam_qc_commands = bismark_bam_qc(output_dir=bam_dir, config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        for cmd in bam_qc_commands:
            command_runner(cmd)

    # summarize bismark bam qc
    summarize_bismark_bam_qc(output_dir=bam_dir)

    # merge bam
    # final_bam_record ['uid', 'index_name', 'bam_path']
    final_bam_record, final_bam_commands = merge_bam(output_dir=bam_dir, record_name='bismark_bam_qc.records.csv')

    # runner
    if qsub:
        raise NotImplementedError
    else:
        for cmd in final_bam_commands:
            command_runner(cmd)

    # generate ALLC file
    generate_allc(output_dir=bam_dir, config=config_path)

    # summarize generate ALLC
    pass
