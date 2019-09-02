import pathlib

from cemba_data.mapping import \
    bismark_mapping, \
    command_runner, \
    summarize_bismark_mapping, \
    bismark_bam_qc, \
    summarize_bismark_bam_qc, \
    merge_bam, \
    generate_allc, summarize_generate_allc, prepare_select_dna_reads, summarize_select_dna_reads


def pipeline_mc(output_dir, config_path, mct=False, qsub=True, cpu=10):
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
        command_runner(bismark_commands, runner=None, cpu=cpu)

    # summarize bismark mapping
    summarize_bismark_mapping(output_dir=bam_dir)

    # prepare bismark bam qc
    # bam_qc_records ['uid', 'index_name', 'read_type', 'bam_path']
    bam_qc_records, bam_qc_commands = bismark_bam_qc(output_dir=bam_dir, config=config_path)

    # runner
    if qsub:
        raise NotImplementedError
    else:
        command_runner(bam_qc_commands, runner=None, cpu=cpu)

    # summarize bismark bam qc
    summarize_bismark_bam_qc(output_dir=bam_dir)

    if mct:
        # prepare select DNA
        select_dna_records, select_dna_commands = prepare_select_dna_reads(output_dir=bam_dir, config=config_path)

        # run select DNA
        if qsub:
            raise NotImplementedError
        else:
            command_runner(select_dna_commands, runner=None, cpu=cpu)

        # summarize select DNA
        summarize_select_dna_reads(output_dir=bam_dir)

        # merge the R1 R2 dna bam
        final_bam_record, final_bam_commands = merge_bam(output_dir=bam_dir,
                                                         record_name='select_dna_reads.records.csv')
    else:
        # merge R1 R2 bam
        # final_bam_record ['uid', 'index_name', 'bam_path']
        final_bam_record, final_bam_commands = merge_bam(output_dir=bam_dir,
                                                         record_name='bismark_bam_qc.records.csv')
    # runner merge bam
    if qsub:
        raise NotImplementedError
    else:
        command_runner(final_bam_commands, runner=None, cpu=cpu)

    # generate ALLC file
    allc_dir = output_dir / 'allc'
    allc_dir.mkdir(exist_ok=True, parents=True)
    allc_record, allc_commands = generate_allc(input_dir=bam_dir,
                                               output_dir=allc_dir,
                                               config=config_path)

    # runner merge bam
    if qsub:
        raise NotImplementedError
    else:
        command_runner(allc_commands, runner=None, cpu=cpu)

    # summarize generate ALLC
    summarize_generate_allc(allc_dir)

    # final summary

