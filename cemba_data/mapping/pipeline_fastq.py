import pathlib

from cemba_data.mapping import \
    make_fastq_dataframe, \
    demultiplex, \
    merge_lane, \
    fastq_qc, \
    command_runner, \
    summarize_demultiplex, \
    summarize_fastq_qc


def pipeline_fastq(input_fastq_pattern, output_dir, config_path, qsub=True):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    fastq_dir.mkdir(exist_ok=True, parents=True)

    # make fastq dataframe
    make_fastq_dataframe(file_path=input_fastq_pattern,
                         output_path=str(fastq_dir / 'fastq_dataframe.csv'),
                         skip_broken_name=False)

    # prepare demultiplex
    # demultiplex_records columns ['uid', 'lane', 'r1_path_pattern', 'r2_path_pattern']
    demultiplex_records, demultiplex_commands = demultiplex(output_dir=fastq_dir, config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        from cemba_data.mapping.demultiplex import demultiplex_runner
        for cmd in demultiplex_commands:
            demultiplex_runner(cmd)

    # summarize demultiplex
    summarize_demultiplex(output_dir=fastq_dir, config=config_path)

    # prepare merge lane
    # merge_lane_records ['uid', 'index_name', 'read_type', 'fastq_path']
    merge_lane_records, merge_lane_commands = merge_lane(output_dir=fastq_dir, config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        for cmd in merge_lane_commands:
            command_runner(cmd)

    # prepare fastq qc
    # fastq_qc_records ['uid', 'index_name', 'read_type', 'fastq_path']
    fastq_qc_records, fastq_qc_commands = fastq_qc(output_dir=fastq_dir, config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        from cemba_data.mapping.fastq_qc import fastq_qc_runner
        for cmd in fastq_qc_commands:
            fastq_qc_runner(cmd)

    # summarize fastq qc
    summarize_fastq_qc(output_dir=fastq_dir)
