import configparser
import logging
import os
import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

from ALLCools._open import open_bam

# logger

log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def get_configuration(config_path=None):
    """
    Read .ini config file from given path
    """
    ref_path_config = configparser.ConfigParser()
    if config_path is None:
        log.info('config path not provided, use default config')
        ref_path_config.read(os.path.dirname(__file__) + '/mapping_config.ini')
    else:
        ref_path_config.read(config_path)
    return ref_path_config


def test_cmd(tool_name, cmd_list):
    try:
        subprocess.run(cmd_list,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       encoding='utf8',
                       check=True)
    except subprocess.CalledProcessError as e:
        log.error(f'Test {tool_name} got non-zero return code {e.returncode}')
        log.error(e.stderr)
        raise
    return


def valid_environments(config):
    log.info('Test mapping environments')

    # test cutadapt
    test_cmd(tool_name='cutadapt', cmd_list=['cutadapt', '--version'])
    # test samtools
    test_cmd(tool_name='samtools', cmd_list=['samtools', '--version'])
    # test picard, picard always have return code 1...
    test_cmd(tool_name='picard', cmd_list=['which', 'picard'])
    # test bismark_mapping
    test_cmd(tool_name='bismark_mapping', cmd_list=['bismark_mapping', '--version'])
    # test bowtie2
    test_cmd(tool_name='bowtie2', cmd_list=['bowtie2', '--version'])
    # test pigz
    test_cmd(tool_name='pigz', cmd_list=['pigz', '-V'])

    bismark_dir = pathlib.Path(config['bismark_mapping']['bismark_reference'])
    if not bismark_dir.is_dir():
        raise TypeError(f"Bismark reference must be a directory contain a sub-dir named Bisulfite_Genome, "
                        f"generated by bismark_genome_preparation. Got a file path")
    if not bismark_dir.exists():
        raise FileNotFoundError(f"Bismark reference directory not found. "
                                f"Path in the config.ini is {bismark_dir}")

    allc_ref_fasta = pathlib.Path(config['callMethylation']['reference_fasta'])
    allc_ref_fai = pathlib.Path(config['callMethylation']['reference_fasta'] + '.fai')
    if not allc_ref_fasta.exists():
        raise FileNotFoundError(f"Reference fasta for ALLC generation not found. "
                                f"Path in the config.ini is {allc_ref_fasta}")
    if not allc_ref_fai.exists():
        raise FileNotFoundError(f".fai index for reference fasta not found. "
                                f"Path of fadix should be {allc_ref_fai}. "
                                f"You can use 'samtools fadix {allc_ref_fasta}' to generate.")
    return


def parse_index_fasta(fasta_path):
    records = {}
    with open(fasta_path) as f:
        key_line = True
        for line in f:
            if key_line:
                key = line.lstrip('>').rstrip('\n')
                key_line = False
            else:
                value = line.lstrip('^').rstrip('\n')
                records[key] = value
                key_line = True
    return records


def command_runner(commands, runner=None, cpu=1):
    if runner is None:
        from functools import partial
        runner = partial(subprocess.run,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         encoding='utf8',
                         shell=True,
                         check=True)
    with ProcessPoolExecutor(cpu) as pool:
        futures = []
        for command in commands:
            future = pool.submit(runner, command)
            futures.append(future)

        for future in as_completed(futures):
            try:
                future.result()
            except subprocess.CalledProcessError as e:
                print("Got error in fastq_qc, command was:")
                print(command)
                print(e.stdout)
                print(e.stderr)
                raise e
    return


def get_bam_header_str(bam_path):
    bam_header = ''
    with open_bam(bam_path, include_header=True) as f:
        for line in f:
            if line.startswith('@'):
                bam_header += line
            else:
                break
    return bam_header
