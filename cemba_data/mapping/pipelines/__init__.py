import os
import pathlib

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


def write_commands(output_dir, cores_per_job, script_dir):
    output_dir_name = output_dir.name
    cmds = {}
    snake_files = list(output_dir.glob('*/Snakefile'))
    for snake_file in snake_files:
        uid = snake_file.parent.name
        cmd = f'snakemake ' \
              f'-d $SCRATCH/{output_dir_name}/{snake_file.parent.name} ' \
              f'--snakefile $SCRATCH/{output_dir_name}/{snake_file.parent.name}/Snakefile ' \
              f'-j {cores_per_job}'
        cmds[uid] = cmd
    uid_order = pd.read_csv(
        output_dir / 'stats/UIDTotalCellInputReadPairs.csv', index_col=0, squeeze=True, header=None
    ).sort_values(ascending=False)
    script_path = script_dir / 'snakemake_cmd.txt'
    with open(script_path, 'w') as f:
        for uid in uid_order.index:
            if uid in cmds:
                f.write(cmds.pop(uid) + '\n')
        try:
            assert len(cmds) == 0
        except AssertionError as e:
            print(cmds)
            print(uid_order)
            raise e
    return f'$SCRATCH/{output_dir_name}/snakemake/sbatch/snakemake_cmd.txt'


def prepare_qsub(name, snakemake_dir, total_jobs, cores_per_job, memory_per_core):
    output_dir = snakemake_dir.parent
    qsub_dir = snakemake_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)
    script_path = write_commands(output_dir, cores_per_job, script_dir=qsub_dir)
    qsub_str = f"""
#!/bin/bash
#$ -N {name}
#$ -V
#$ -l h_rt=999:99:99
#$ -l s_rt=999:99:99
#$ -wd {qsub_dir}
#$ -e {qsub_dir}/qsub.error.log
#$ -o {qsub_dir}/qsub.output.log
#$ -pe smp 1
#$ -l h_vmem=3G

yap qsub \
--command_file_path {script_path} \
--working_dir {qsub_dir} \
--project_name {name} \
--total_cpu {int(cores_per_job * total_jobs)} \
--qsub_global_parms "-pe smp={cores_per_job};-l h_vmem={memory_per_core}"
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


def prepare_sbatch(name, snakemake_dir, sbatch_cores_per_job=38):
    output_dir = snakemake_dir.parent
    output_dir_name = output_dir.name
    mode = get_configuration(output_dir / 'mapping_config.ini')['mode']
    if mode == 'm3c':
        print('Decrease cores_per_job to 29 for m3C library.')
        sbatch_cores_per_job = 29
    sbatch_dir = snakemake_dir / 'sbatch'
    sbatch_dir.mkdir(exist_ok=True)

    script_path = write_commands(output_dir, cores_per_job=sbatch_cores_per_job, script_dir=sbatch_dir)
    # the path here is using stampede path
    sbatch_cmd = f'yap sbatch ' \
                 f'--project_name {name} ' \
                 f'--command_file_path {script_path} ' \
                 f'--working_dir $SCRATCH/{output_dir_name}/snakemake/sbatch ' \
                 f'--time_str 12:00:00'  # TODO better estimate time based on total reads
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
          f"Alternatively, you can sbatch the commands in {script_path} by yourself, "
          f"as long as they all get successfully executed.")
    print('#' * 60 + '\n')
    return


def prepare_run(output_dir, total_jobs=12, cores_per_job=8, memory_per_core='5G', name=None):
    # TODO test what parameter works best in real nova-seq data
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
                     memory_per_core=memory_per_core)
        prepare_sbatch(name=name, snakemake_dir=snakemake_dir)
    else:
        script_path = write_commands(output_dir, cores_per_job, script_dir=snakemake_dir)
        print(f"All snakemake commands need to be executed were summarized in {script_path}")
        print(f"You need to execute them based on the computational environment you have "
              f"(e.g., use a job scheduler or run locally).")

    print(f"Once all commands are executed successfully, use 'yap summary' to generate final mapping summary.")
    return
