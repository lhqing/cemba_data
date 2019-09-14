import pathlib

from cemba_data.mapping import \
    bismark_mapping, \
    command_runner, \
    bismark_bam_qc, \
    merge_bam, \
    generate_allc, prepare_select_dna_reads
from cemba_data.qsub import qsub


def pipeline_mc(output_dir, config_path, mct=False, mode='command_only', cpu=10):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'bismark_bam'
    bam_dir.mkdir(exist_ok=True, parents=True)
    allc_dir = output_dir / 'allc'
    allc_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub/bismark_bam_allc'
    qsub_dir.mkdir(exist_ok=True, parents=True)

    # prepare bismark mapping
    # bismark_records ['uid', 'index_name', 'read_type', 'bam_path']
    bismark_records, bismark_commands = bismark_mapping(input_dir=fastq_dir,
                                                        output_dir=bam_dir,
                                                        config=config_path)
    with open(qsub_dir / 'bismark_commands.txt', 'w') as f:
        f.write('\n'.join(bismark_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'bismark_commands.txt'),
             working_dir=qsub_dir,
             project_name='bismark',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp 4;-l h_vmem=5G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(bismark_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # prepare bismark bam qc
    # bam_qc_records ['uid', 'index_name', 'read_type', 'bam_path']
    bam_qc_records, bam_qc_commands = bismark_bam_qc(output_dir=bam_dir, config=config_path)
    with open(qsub_dir / 'bam_qc_commands.txt', 'w') as f:
        f.write('\n'.join(bam_qc_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'bam_qc_commands.txt'),
             working_dir=qsub_dir,
             project_name='bam_qc',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp 2;-l h_vmem=4G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(bam_qc_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    if mct:
        # prepare select DNA
        select_dna_records, select_dna_commands = prepare_select_dna_reads(output_dir=bam_dir, config=config_path)
        with open(qsub_dir / 'select_dna_commands.txt', 'w') as f:
            f.write('\n'.join(select_dna_commands))

        # runner
        if mode == 'qsub':
            qsub(command_file_path=str(qsub_dir / 'select_dna_commands.txt'),
                 working_dir=qsub_dir,
                 project_name='select_dna',
                 wait_until=None,
                 total_cpu=cpu,
                 total_mem=500,
                 force_redo=False,
                 qsub_global_parms='-pe smp 1;-l h_vmem=5G')
        elif mode == 'command_only':
            pass
        elif mode == 'local':
            command_runner(select_dna_commands, runner=None, cpu=cpu)
        else:
            raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')
        records_file = 'select_dna_reads.records.csv'
    else:
        records_file = 'bismark_bam_qc.records.csv'

    # merge R1 R2 bam
    # final_bam_record ['uid', 'index_name', 'bam_path']
    final_bam_record, final_bam_commands = merge_bam(output_dir=bam_dir,
                                                     record_name=records_file)
    with open(qsub_dir / 'final_bam_commands.txt', 'w') as f:
        f.write('\n'.join(final_bam_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'final_bam_commands.txt'),
             working_dir=qsub_dir,
             project_name='merge_bam',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp 1;-l h_vmem=5G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(final_bam_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # generate ALLC file
    allc_record, allc_commands = generate_allc(input_dir=bam_dir,
                                               output_dir=allc_dir,
                                               config=config_path)
    with open(qsub_dir / 'allc_commands.txt', 'w') as f:
        f.write('\n'.join(allc_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'allc_commands.txt'),
             working_dir=qsub_dir,
             project_name='allc',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp 2;-l h_vmem=4G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(allc_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')
