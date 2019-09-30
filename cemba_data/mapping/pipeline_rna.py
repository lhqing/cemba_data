import json
import pathlib

from cemba_data.mapping import \
    prepare_select_rna_reads, \
    command_runner, \
    star_mapping, \
    star_bam_qc, \
    prepare_feature_count
from cemba_data.qsub import qsub


def pipeline_rna(output_dir, config_path, mode='command_only', cpu=10):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'star_bam'
    bam_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub/star_bam'
    qsub_dir.mkdir(exist_ok=True, parents=True)

    # prepare STAR mapping
    star_mapping_records, star_mapping_commands = star_mapping(input_dir=fastq_dir,
                                                               output_dir=bam_dir,
                                                               config=config_path)
    with open(qsub_dir / 'star_mapping_commands.txt', 'w') as f:
        f.write('\n'.join(star_mapping_commands))
    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'star_mapping_commands.txt'),
             working_dir=qsub_dir,
             project_name='star_mapping',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp=6;-l h_vmem=6G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(star_mapping_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # prepare STAR BAM QC
    star_bam_qc_records, star_bam_qc_commands = star_bam_qc(output_dir=bam_dir, config=config_path)
    with open(qsub_dir / 'star_bam_qc_commands.txt', 'w') as f:
        f.write('\n'.join(star_bam_qc_commands))
    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'star_bam_qc_commands.txt'),
             working_dir=qsub_dir,
             project_name='star_bam_qc',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp=1;-l h_vmem=4G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(star_bam_qc_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # prepare select RNA
    select_rna_records, select_rna_commands = prepare_select_rna_reads(output_dir=bam_dir, config=config_path)
    with open(qsub_dir / 'select_rna_commands.txt', 'w') as f:
        f.write('\n'.join(select_rna_commands))
    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'select_rna_commands.txt'),
             working_dir=qsub_dir,
             project_name='select_rna',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp=1;-l h_vmem=4G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(select_rna_commands, runner=None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # feature count
    feature_count_output, feature_count_command = prepare_feature_count(output_dir=bam_dir, config=config_path)
    with open(qsub_dir / 'feature_count_commands.txt', 'w') as f:
        f.write(feature_count_command)
    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'feature_count_commands.txt'),
             working_dir=qsub_dir,
             project_name='feature_count',
             wait_until=None,
             total_cpu=10,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp=10;-l h_vmem=5G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner([feature_count_command], runner=None, cpu=1)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')
