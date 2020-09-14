import os
import pathlib

import cemba_data
from .mc import mc_config_str
from ...utilities import get_configuration

# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


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


def prepare_run(output_dir, total_jobs, cores_per_job, memory_per_core, name=None):
    if cores_per_job < 4:
        raise ValueError(f'cores must >= 4 to run this pipeline.')
    output_dir = pathlib.Path(output_dir).absolute()
    if name is None:
        name = output_dir.name
    snakemake_dir = output_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True)

    cmds = []
    snake_files = list(output_dir.glob('*/Snakefile'))
    for snake_file in snake_files:
        cmd = f'snakemake -d {snake_file.parent} --snakefile {snake_file} -j {cores_per_job}'
        cmds.append(cmd)
    with open(snakemake_dir / 'snakemake_cmd.txt', 'w') as f:
        f.write('\n'.join(cmds))

    # this is only some automatic code for ecker lab...
    # so conditioned by the host name
    host_name = os.environ['HOSTNAME']
    if host_name[:4] in ['bpho', 'gale']:
        qsub_str = f"""
        #!/bin/bash
        #$ -N {name}
        #$ -V
        #$ -l h_rt=999:99:99
        #$ -l s_rt=999:99:99
        #$ -wd {snakemake_dir}
        #$ -e {snakemake_dir}/qsub.error.log
        #$ -o {snakemake_dir}/qsub.output.log
        #$ -pe smp 1
        #$ -l h_vmem=3G

        yap qsub \
        --command_file_path {snakemake_dir}/snakemake_cmd.txt \
        --working_dir {snakemake_dir} \
        --project_name {name} \
        --total_cpu {int(cores_per_job * total_jobs)} \
        --qsub_global_parms "-pe smp={cores_per_job};-l h_vmem={memory_per_core}"
        """

        with open(snakemake_dir / 'qsub.sh', 'w') as f:
            f.write(qsub_str)
        print(f"All snakemake commands need to be executed were summarized in {snakemake_dir / 'qsub.sh'}")
        print(f'You just need to qsub this script to map the whole library in {output_dir}')
    else:
        print(f"All snakemake commands need to be executed were summarized in {snakemake_dir / 'snakemake_cmd.txt'}")
        print(f"You need to execute them based on the computational environment you have "
              f"(e.g., use a job scheduler or run locally).")

    print(f"Once all commands are executed successfully, use 'yap summary' to generate final mapping summary.")
    return
