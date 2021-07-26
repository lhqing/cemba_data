import pathlib
import pandas as pd
from .prepare_impute import execute_command


def prepare_dataset_commands(output_dir, fasta_path, cpu=10):
    output_dir = pathlib.Path(output_dir)
    project_name = output_dir.name
    scool_dir = output_dir / 'scool'
    snakemake_dir = scool_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True, parents=True)
    raw_dir = scool_dir / 'raw'
    raw_dir.mkdir(exist_ok=True)
    impute_dir = scool_dir / 'impute'
    impute_dir.mkdir(exist_ok=True)
    dataset_dir = scool_dir / 'dataset'
    dataset_dir.mkdir(exist_ok=True)

    # Calculate compartment at 100Kb resolution
    compartment_input_dir = impute_dir / '100K'
    compartment_cell_table = pd.Series({
        path.name.split('.')[0]: str(path)
        for path in compartment_input_dir.glob('*/*.cool')
    })
    compartment_cell_table_path = compartment_input_dir / 'cell_table.tsv'
    compartment_cell_table.to_csv(compartment_cell_table_path, sep='\t', header=None)
    # prepare a whole genome CpG ratio profile
    cpg_path = compartment_input_dir / 'cpg_ratio.hdf'
    cpg_ratio_cmd = f'hicluster cpg-ratio --cell_url {compartment_cell_table.iloc[0]} ' \
                    f'--fasta_path {fasta_path} --hdf_output_path {cpg_path}'
    execute_command(cpg_ratio_cmd)
    # compartment command
    compartment_cmd = f'hicluster compartment ' \
                      f'--cell_table_path {compartment_cell_table_path} ' \
                      f'--output_prefix {dataset_dir / project_name} ' \
                      f'--cpg_profile_path {cpg_path} ' \
                      f'--cpu {cpu}'

    # Calculate domain at 25Kb resolution
    domain_input_dir = impute_dir / '25K'
    domain_cell_table = pd.Series({
        path.name.split('.')[0]: str(path)
        for path in domain_input_dir.glob('*/*.cool')
    })
    domain_cell_table_path = domain_input_dir / 'cell_table.tsv'
    domain_cell_table.to_csv(domain_cell_table_path, sep='\t', header=None)
    domain_cmd = f'hicluster domain ' \
                 f'--cell_table_path {domain_cell_table_path} ' \
                 f'--output_prefix {dataset_dir / project_name} ' \
                 f'--resolution 25000 ' \
                 f'--window_size 10 ' \
                 f'--cpu {cpu}'

    # Calculate cell embedding/decomposition at 100Kb resolution
    embedding_dir = dataset_dir / 'embedding'
    embedding_dir.mkdir(exist_ok=True)
    embedding_cmd = f'hicluster embedding ' \
                    f'--cell_table_path {compartment_cell_table_path} ' \
                    f'--output_dir {embedding_dir} ' \
                    f'--dim 50 ' \
                    f'--dist 1000000 ' \
                    f'--resolution 100000 ' \
                    f'--scale_factor 100000 ' \
                    f'--norm_sig --save_raw ' \
                    f'--cpu {cpu}'

    # prepare qsub
    qsub_dir = snakemake_dir / 'qsub'
    with open(qsub_dir / 'dataset_cmd.txt', 'w') as f:
        f.write('\n'.join([compartment_cmd, domain_cmd, embedding_cmd]))
    qsub_str = f"""
#!/bin/bash
#$ -N y{project_name}
#$ -V
#$ -l h_rt=999:99:99
#$ -l s_rt=999:99:99
#$ -wd {qsub_dir}
#$ -e {qsub_dir}/qsub_dataset.error.log
#$ -o {qsub_dir}/qsub_dataset.output.log
#$ -pe smp 1
#$ -l h_vmem=3G

yap qsub --command_file_path {qsub_dir}/dataset_cmd.txt \
--working_dir {qsub_dir} --project_name y{project_name}_dataset \
--total_cpu {int(cpu*3)} --qsub_global_parms "-pe smp={cpu};-l h_vmem=5G"
"""
    with open(qsub_dir / 'qsub_dataset.sh', 'w') as f:
        f.write(qsub_str)
    return
