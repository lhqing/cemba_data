import glob
import os
import pathlib
import subprocess

import pandas as pd

import cemba_data
from .m3c import m3c_config_str
from .mc import mc_config_str
from .mct import mct_config_str
from ...utilities import get_configuration

# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])
INHOUSE_SERVERS = ['bpho', 'gale', 'cemba', 'oberon']


def prepare_uid_snakefile(uid_dir, config_str, snake_template):
    cell_ids = [path.name.split('.')[0][:-3] for path in (uid_dir / 'fastq').glob('*R1.fq.gz')]
    cell_id_str = f'CELL_IDS = {cell_ids}\n'

    # no file in this UID, do not make snakefile
    if len(cell_ids) == 0:
        return

    total_snakefile = config_str + cell_id_str + snake_template
    with open(uid_dir / 'Snakefile', 'w') as f:
        f.write(total_snakefile)
    return


def validate_mapping_config(output_dir):
    output_dir = pathlib.Path(output_dir).absolute()
    config = get_configuration(output_dir / 'mapping_config.ini')
    try:
        mode = config['mode']
    except KeyError:
        raise KeyError('mode not found in the config file.')

    if mode == 'mc':
        config_str = mc_config_str(config)
    elif mode == 'mct':
        config_str = mct_config_str(config)
    elif mode == 'm3c':
        config_str = m3c_config_str(config)
    else:
        raise ValueError(f'Unknown mode {mode}')

    print(f'Mapping config file looks good. Here is what will be used in generating Snakefile:\n{config_str}')
    return


def make_snakefile(output_dir):
    output_dir = pathlib.Path(output_dir).absolute()
    config = get_configuration(output_dir / 'mapping_config.ini')
    try:
        mode = config['mode']
    except KeyError:
        raise KeyError('mode not found in the config file.')

    if mode == 'mc':
        config_str = mc_config_str(config)
    elif mode == 'mct':
        config_str = mct_config_str(config)
    elif mode == 'm3c':
        config_str = m3c_config_str(config)
    else:
        raise ValueError(f'Unknown mode {mode}')
    print('Making Snakefile based on mapping config INI file. The parameters are:')
    print(config_str)

    with open(PACKAGE_DIR / f'mapping/Snakefile_template/{mode}.Snakefile') as f:
        snake_template = f.read()

    for sub_dir in output_dir.iterdir():
        if sub_dir.is_dir():
            if sub_dir.name not in ['stats', 'snakemake']:
                prepare_uid_snakefile(uid_dir=sub_dir,
                                      config_str=config_str,
                                      snake_template=snake_template)
    return


def write_qsub_commands(output_dir, cores_per_job, memory_gb_per_core, script_dir):
    memory_per_core = int(memory_gb_per_core[:-1]) * 1000
    cmds = {}
    snake_files = list(output_dir.glob('*/Snakefile'))
    for snake_file in snake_files:
        uid = snake_file.parent.name
        cmd = f'snakemake ' \
              f'-d {snake_file.parent} ' \
              f'--snakefile {snake_file} ' \
              f'-j {cores_per_job} ' \
              f'--default-resources mem_mb=100 ' \
              f'--resources mem_mb={int(cores_per_job * memory_per_core)} '
        cmds[uid] = cmd
    script_path = script_dir / 'snakemake_cmd.txt'
    with open(script_path, 'w') as f:
        try:
            uid_order = pd.read_csv(
                output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0, squeeze=True, header=None
            ).sort_values(ascending=False)
            for uid in uid_order.index:
                if uid in cmds:
                    f.write(cmds.pop(uid) + '\n')
            try:
                assert len(cmds) == 0
            except AssertionError as e:
                print(cmds)
                print(uid_order)
                raise e
        except FileNotFoundError:
            # uid_order file do not exist (when starting from cell FASTQs)
            for cmd in cmds.values():
                f.write(cmd + '\n')
    return script_path


def write_sbatch_commands(output_dir, cores_per_job, script_dir, total_mem_mb):
    output_dir_name = output_dir.name
    cmds = {}
    snake_files = list(output_dir.glob('*/Snakefile'))
    for snake_file in snake_files:
        uid = snake_file.parent.name
        cmd = f'snakemake ' \
              f'-d $SCRATCH/{output_dir_name}/{snake_file.parent.name} ' \
              f'--snakefile $SCRATCH/{output_dir_name}/{snake_file.parent.name}/Snakefile ' \
              f'-j {cores_per_job} ' \
              f'--default-resources mem_mb=100 ' \
              f'--resources mem_mb={total_mem_mb} '
        cmds[uid] = cmd
    script_path = script_dir / 'snakemake_cmd.txt'
    with open(script_path, 'w') as f:
        try:
            uid_order = pd.read_csv(
                output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0, squeeze=True, header=None
            ).sort_values(ascending=False)
            for uid in uid_order.index:
                if uid in cmds:
                    f.write(cmds.pop(uid) + '\n')
            try:
                assert len(cmds) == 0
            except AssertionError as e:
                print(cmds)
                print(uid_order)
                raise e
        except FileNotFoundError:
            # uid_order file do not exist (when starting from cell FASTQs)
            for cmd in cmds.values():
                f.write(cmd + '\n')
    return f'$SCRATCH/{output_dir_name}/snakemake/sbatch/snakemake_cmd.txt'


def prepare_qsub(name, snakemake_dir, total_jobs, cores_per_job, memory_gb_per_core):
    output_dir = snakemake_dir.parent
    qsub_dir = snakemake_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)
    script_path = write_qsub_commands(output_dir, cores_per_job, memory_gb_per_core, script_dir=qsub_dir)
    qsub_str = f"""
#!/bin/bash
#$ -N yap{name}
#$ -V
#$ -l h_rt=99:99:99
#$ -l s_rt=99:99:99
#$ -wd {qsub_dir}
#$ -e {qsub_dir}/qsub.error.log
#$ -o {qsub_dir}/qsub.output.log
#$ -pe smp 1
#$ -l h_vmem=3G

yap qsub \
--command_file_path {script_path} \
--working_dir {qsub_dir} \
--project_name y{name} \
--total_cpu {int(cores_per_job * total_jobs)} \
--qsub_global_parms "-pe smp={cores_per_job};-l h_vmem={memory_gb_per_core}"
"""
    qsub_total_path = qsub_dir / 'qsub.sh'
    with open(qsub_total_path, 'w') as f:
        f.write(qsub_str)
    print('#' * 60)
    print(f"IF YOU USE QSUB ON GALE: ")
    print(f"All snakemake commands need to be executed "
          f"were included in {qsub_total_path}")
    print(f"You just need to qsub this script to "
          f"map the whole library in {output_dir}")
    print(f"You can also change the per job parameters in {script_path} "
          f"or change the global parameters in {qsub_total_path}")
    print(f"Read 'yap qsub -h' if you want to have more options about sbatch. "
          f"Alternatively, you can sbatch the commands in {script_path} by yourself, "
          f"as long as they all get successfully executed.")
    print('#' * 60 + '\n')
    return


def prepare_sbatch(name, snakemake_dir, sbatch_cores_per_job=96, total_mem_mb=192000):
    output_dir = snakemake_dir.parent
    output_dir_name = output_dir.name
    mode = get_configuration(output_dir / 'mapping_config.ini')['mode']
    time_str = "8:00:00"
    if mode == 'm3c':
        time_str = "10:00:00"
    elif mode == 'mc':
        time_str = "7:00:00"
    elif mode == 'mct':
        time_str = "8:00:00"
    sbatch_dir = snakemake_dir / 'sbatch'
    sbatch_dir.mkdir(exist_ok=True)

    script_path = write_sbatch_commands(output_dir,
                                        cores_per_job=sbatch_cores_per_job,
                                        script_dir=sbatch_dir,
                                        total_mem_mb=total_mem_mb)
    # the path here is using stampede path
    sbatch_cmd = f'yap sbatch ' \
                 f'--project_name {name} ' \
                 f'--command_file_path {script_path} ' \
                 f'--working_dir $SCRATCH/{output_dir_name}/snakemake/sbatch ' \
                 f'--time_str {time_str}'
    sbatch_total_path = sbatch_dir / 'sbatch.sh'
    with open(sbatch_total_path, 'w') as f:
        f.write(sbatch_cmd)
    print('#' * 60)
    print(f"IF YOU USE SBATCH ON STAMPEDE2:")
    print(f"All snakemake commands need to be executed "
          f"were included in {sbatch_total_path}")
    print(f"You just need to run this script to "
          f"map the whole library in {output_dir}. "
          f"Note that this script will not return until all the mapping job finished. "
          f"So you better run this script with nohup or in a screen.")
    print(f"You can also change "
          f"the per job parameters in {script_path} "
          f"or change the global parameters in {sbatch_total_path}")
    print(f"Read 'yap sbatch -h' if you want to have more options about sbatch. "
          f"Alternatively, you can sbatch the commands in "
          f"$SCRATCH/{output_dir_name}/snakemake/sbatch/sbatch.sh by yourself, "
          f"as long as they all get successfully executed.")
    print('#' * 60 + '\n')
    return


def prepare_run(output_dir, total_jobs=12, cores_per_job=10, memory_gb_per_core='5G', name=None):
    config = get_configuration(output_dir / 'mapping_config.ini')
    mode = config['mode']
    if mode in ['mc', 'm3c'] and cores_per_job < 4:
        raise ValueError(f'cores must >= 4 to run this pipeline.')
    elif mode == 'mct' and cores_per_job < 10:
        raise ValueError(f'cores must >= 10 to run this pipeline.')

    output_dir = pathlib.Path(output_dir).absolute()
    if name is None:
        name = output_dir.name
    snakemake_dir = output_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True)

    # this is only some automatic code for ecker lab...
    # so conditioned by the host name
    host_name = os.environ['HOSTNAME']
    if any([host_name.startswith(s) for s in INHOUSE_SERVERS]):
        prepare_qsub(name=name,
                     snakemake_dir=snakemake_dir,
                     total_jobs=total_jobs,
                     cores_per_job=cores_per_job,
                     memory_gb_per_core=memory_gb_per_core)
        prepare_sbatch(name=name, snakemake_dir=snakemake_dir)
    else:
        script_path = write_qsub_commands(output_dir, cores_per_job, memory_gb_per_core, script_dir=snakemake_dir)
        print(f"All snakemake commands need to be executed were summarized in {script_path}")
        print(f"You need to execute them based on the computational environment you have "
              f"(e.g., use a job scheduler or run locally).")

    print(f"Once all commands are executed successfully, use 'yap summary' to generate final mapping summary.")
    return


def start_from_cell_fastq(output_dir, fastq_pattern, config_path):
    output_dir = pathlib.Path(output_dir).absolute()
    if output_dir.exists():
        raise FileExistsError(f'Output dir {output_dir} already exist, please delete it or use another path.')
    output_dir.mkdir()
    subprocess.run(['cp', config_path, f'{output_dir}/mapping_config.ini'], check=True)
    stats_dir = output_dir / 'stats'
    stats_dir.mkdir(exist_ok=True)

    # parse fastq patterns
    fastq_paths = [pathlib.Path(p).absolute() for p in glob.glob(fastq_pattern)]
    r1_records = {}
    r2_records = {}
    for path in fastq_paths:
        *cell_id, suffix = path.name.split('-')
        cell_id = '-'.join(cell_id)
        if suffix == 'R1.fq.gz':
            if cell_id in r1_records:
                raise ValueError(f'Found duplicated cell ID: {cell_id}, '
                                 f'File caused this error: {path}')
            r1_records[cell_id] = path
        elif suffix == 'R2.fq.gz':
            if cell_id in r2_records:
                raise ValueError(f'Found duplicated cell ID: {cell_id}, '
                                 f'File caused this error: {path}')
            r2_records[cell_id] = path
        else:
            raise ValueError(
                f'Unable to parse read type. Expect file name ends with "-R1.fq.gz" or "-R2.fq.gz", '
                f'File caused this error: {path}'
            )
    fastq_df = pd.DataFrame({'R1Path': r1_records, 'R2Path': r2_records})

    # make symlink of fastq files, using dir structure of demultiplex
    # the cells are randomly grouped though, max group is 64
    groups = min(64, fastq_df.shape[0])
    for i, (cell_id, (r1_path, r2_path)) in enumerate(fastq_df.sample(fastq_df.shape[0]).iterrows()):
        group_id = i % groups
        fastq_dir = output_dir / f'Group{group_id}/fastq'
        fastq_dir.mkdir(exist_ok=True, parents=True)

        # make symlinks
        new_r1_path = fastq_dir / r1_path.name
        new_r1_path.symlink_to(r1_path)
        new_r2_path = fastq_dir / r2_path.name
        new_r2_path.symlink_to(r2_path)

    # prepare scripts
    make_snakefile(output_dir)
    prepare_run(output_dir)
    return
