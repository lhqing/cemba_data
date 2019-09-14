import logging
import os
import pathlib
import subprocess

from cemba_data.mapping import get_configuration

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

MAPPING_MODE_CHOICES = ['mc', 'mct', 'nome', 'mct-nome']


def print_default_configuration(mode='mc'):
    """
    Print default .ini config file or save to out_path
    """
    mode = mode.lower()
    if mode not in MAPPING_MODE_CHOICES:
        raise ValueError(f'Unknown mode: {mode}')
    with open(os.path.dirname(__file__) + f'/mapping_config_{mode}.ini') as f:
        configs = f.readlines()
        for line in configs:
            print(line, end='')
    return


def command_order(output_dir):
    command_files = {path.name.split('.')[0]: path
                     for path in pathlib.Path(output_dir).glob('qsub/*/*command*')}
    command_order_full = ['demultiplex_commands',
                          'merge_lane_commands',
                          'fastq_qc_commands',
                          'bismark_commands',
                          'bam_qc_commands',
                          'select_dna_commands',
                          'final_bam_commands',
                          'star_mapping_commands',
                          'star_bam_qc_commands',
                          'select_rna_commands']
    with open(output_dir / 'command_order.txt', 'w') as f:
        for command_name in command_order_full:
            try:
                f.write(str(command_files[command_name]) + '\n')
            except KeyError:
                continue
    return


def pipeline(input_fastq_pattern,
             output_dir,
             config_path,
             fastq_dataframe_path=None,
             mode='command_only',
             cpu=10):
    _output_dir = pathlib.Path(output_dir)
    _output_dir.mkdir(exist_ok=True, parents=True)
    _config_path = str(_output_dir / pathlib.Path(config_path).name)
    subprocess.run(['cp', str(config_path), _config_path], check=True)
    config_path = _config_path

    if mode == 'command_only':
        if cpu > 10:
            print('Local mode is only for testing, better set a smaller cpu. Change to 10.')
            cpu = 10

    config = get_configuration(config_path)
    mct = 'mct' in config['mode']['mode'].lower()

    # test environment
    from cemba_data.mapping.test_environment import testing_mapping_installation
    testing_mapping_installation(mct=mct)

    # pipeline_fastq
    from cemba_data.mapping.pipeline_fastq import pipeline_fastq
    pipeline_fastq(input_fastq_pattern,
                   output_dir,
                   config_path,
                   fastq_dataframe_path=fastq_dataframe_path,
                   mode=mode,
                   cpu=cpu)
    # TODO add pipeline that start from post demultiplex

    # pipeline_mc
    from cemba_data.mapping.pipeline_mc import pipeline_mc
    pipeline_mc(output_dir,
                config_path,
                mct=mct,
                mode=mode,
                cpu=cpu)

    # if mct: pipeline_rna
    if mct:
        from cemba_data.mapping.pipeline_rna import pipeline_rna
        pipeline_rna(output_dir,
                     config_path,
                     mode=mode,
                     cpu=cpu)

    if mode == 'command_only':
        command_order(output_dir)
    # TODO each individual pipeline step should be able to run separately
    return 0
