import pathlib
import subprocess

from cemba_data.mapping import \
    make_fastq_dataframe, \
    demultiplex, \
    merge_lane, \
    fastq_qc, \
    command_runner
from cemba_data.qsub import qsub


def pipeline_fastq(input_fastq_pattern,
                   output_dir,
                   config_path,
                   fastq_dataframe_path=None,
                   mode='command_only',
                   cpu=10):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    fastq_dir.mkdir(exist_ok=True, parents=True)
    qsub_dir = output_dir / 'qsub/fastq'
    qsub_dir.mkdir(exist_ok=True, parents=True)

    # STEP 1: make fastq dataframe
    if fastq_dataframe_path is None:
        make_fastq_dataframe(file_path=input_fastq_pattern,
                             output_path=str(fastq_dir / 'fastq_dataframe.csv'),
                             config_path=config_path)
    else:
        subprocess.run(['cp', fastq_dataframe_path, str(fastq_dir / 'fastq_dataframe.csv')])

    # STEP 2: demultiplex
    # prepare demultiplex
    # demultiplex_records columns ['uid', 'lane', 'r1_path_pattern', 'r2_path_pattern']
    demultiplex_records, demultiplex_commands = demultiplex(output_dir=fastq_dir, config=config_path)
    with open(qsub_dir / 'demultiplex_commands.txt', 'w') as f:
        f.write('\n'.join(demultiplex_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'demultiplex_commands.txt'),
             working_dir=qsub_dir,
             project_name='demultiplex',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms='-pe smp=2;-l h_vmem=3G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(demultiplex_commands, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # STEP 3: merge lane
    # prepare merge lane
    # merge_lane_records ['uid', 'index_name', 'read_type', 'fastq_path']
    merge_lane_records, merge_lane_commands = merge_lane(output_dir=fastq_dir)
    with open(qsub_dir / 'merge_lane_commands.txt', 'w') as f:
        f.write('\n'.join(merge_lane_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'merge_lane_commands.txt'),
             working_dir=qsub_dir,
             project_name='merge_lane',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms=f'-pe smp=2;-l h_vmem=3G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(merge_lane_commands, None, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')

    # STEP 4: fastq qc
    # prepare fastq qc
    # fastq_qc_records ['uid', 'index_name', 'read_type', 'fastq_path']
    fastq_qc_records, fastq_qc_commands = fastq_qc(output_dir=fastq_dir, config=config_path)
    with open(qsub_dir / 'fastq_qc_commands.txt', 'w') as f:
        f.write('\n'.join(fastq_qc_commands))

    # runner
    if mode == 'qsub':
        qsub(command_file_path=str(qsub_dir / 'fastq_qc_commands.txt'),
             working_dir=qsub_dir,
             project_name='fastq_qc',
             wait_until=None,
             total_cpu=cpu,
             total_mem=500,
             force_redo=False,
             qsub_global_parms=f'-pe smp=2;-l h_vmem=3G')
    elif mode == 'command_only':
        pass
    elif mode == 'local':
        command_runner(fastq_qc_commands, cpu=cpu)
    else:
        raise ValueError(f'mode can only be in ["qsub", "command_only", "local"], got {mode}')
    return
