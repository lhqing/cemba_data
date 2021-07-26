import pathlib
import pandas as pd
import subprocess
import re


def execute_command(cmd):
    print(cmd)
    try:
        subprocess.run(cmd,
                       stderr=subprocess.PIPE,
                       encoding='utf8',
                       stdout=subprocess.PIPE,
                       check=True,
                       shell=True)
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e
    return


def modify_snakefile_for_sbatch(path, project_name):
    """Change the input_scool line in snakefile to remote scratch path"""
    with open(path) as f:
        lines = f.readlines()

    with open(path, 'w') as f:
        for line in lines:
            if line.startswith('input_scool'):
                original_path = line.strip().strip("'").split(' = ')[1]
                original_path = pathlib.Path(original_path)
                new_path = f"f'{{scratch_path}}/{project_name}/scool/raw/{original_path.name}'"
                new_line = f'''
import subprocess
scratch_path = subprocess.run('echo $SCRATCH', 
                              shell=True,  
                              check=True,  
                              stdout=subprocess.PIPE,  
                              encoding='utf8').stdout.strip()
input_scool = {new_path}
'''
                f.write(new_line)
            else:
                f.write(line)
    return


def prepare_impute_dir(output_dir,
                       chrom_size_path,
                       scheduler='sbatch',
                       scool_cpu=5,
                       cpu_per_job=10,
                       batch_size=100,
                       total_cpu=100,
                       min_contacts_per_cell=50000,
                       sbatch_time_str='2:00:00',
                       min_cutoff='1e-5',
                       skip_scool_prep=False):
    # set up directory structure under output_dir
    output_dir = pathlib.Path(output_dir).absolute()
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

    # make contacts
    contact_paths = pd.Series({path.name.split('.')[0]: str(path)
                               for path in output_dir.glob('*/hic/*contact.tsv.gz')})

    # filter by minimum contacts cutoff
    cell_meta = pd.read_csv(f'{output_dir}/stats/MappingSummary.csv.gz', index_col=0)
    use_cells = cell_meta[cell_meta['CisLongContact'] > min_contacts_per_cell].index
    print(f'{cell_meta.shape[0]} cells in total, '
          f'{use_cells.size} cells having contacts > {min_contacts_per_cell} cutoff.')
    contact_paths = contact_paths.loc[use_cells].copy()

    contact_table_path = scool_dir / 'contacts_table.tsv'
    contact_paths.to_csv(contact_table_path, sep='\t', header=None)

    # target scool file paths
    scool_path_10k = raw_dir / f'{project_name}.10K.scool'
    scool_path_25k = raw_dir / f'{project_name}.25K.scool'
    scool_path_100k = raw_dir / f'{project_name}.100K.scool'

    # generate scool command
    if skip_scool_prep:
        if not scool_path_10k.exists():
            raise FileNotFoundError(f'skip_scool_prep is True, but {scool_path_10k} not found.')
        if not scool_path_25k.exists():
            raise FileNotFoundError(f'skip_scool_prep is True, but {scool_path_25k} not found.')
        if not scool_path_100k.exists():
            raise FileNotFoundError(f'skip_scool_prep is True, but {scool_path_100k} not found.')
    else:
        scool_cmd = f'hicluster generate-scool --contacts_table {contact_table_path} ' \
                f'--output_prefix {raw_dir / project_name} ' \
                f'--chrom_size_path {chrom_size_path} ' \
                f'--resolutions 10000 25000 100000 ' \
                f'--min_pos_dist 2500 --cpu {scool_cpu}'
        print('Generating raw scool files for 10Kb, 25Kb and 100Kb resolutions.')
        execute_command(scool_cmd)

    # prepare imputation commands
    # Prepare 10K impute
    this_out_dir = pathlib.Path(f'{impute_dir}/10K/')
    this_out_dir.mkdir(exist_ok=True)
    impute_10k_cmd = f'hicluster imputation ' \
                     f'--input_scool {scool_path_10k} ' \
                     f'--output_dir {this_out_dir} ' \
                     f'--chrom_size_path {chrom_size_path} ' \
                     f'--output_dist 5050000 ' \
                     f'--window_size 30000000 ' \
                     f'--step_size 10000000 ' \
                     f'--resolution 10000 ' \
                     f'--batch_size {batch_size} ' \
                     f'--pad 2 --std 1 --rp 0.5 --tol 0.01 ' \
                     f'--min_cutoff {min_cutoff} ' \
                     f'--cpu_per_job {cpu_per_job}'
    print('Prepare imputation snakefiles for 10Kb resolution.')
    execute_command(impute_10k_cmd)

    # 25K
    this_out_dir = pathlib.Path(f'{impute_dir}/25K/')
    this_out_dir.mkdir(exist_ok=True)
    impute_25k_cmd = f'hicluster imputation ' \
                     f'--input_scool {scool_path_25k} ' \
                     f'--output_dir {this_out_dir} ' \
                     f'--chrom_size_path {chrom_size_path} ' \
                     f'--output_dist 10050000 ' \
                     f'--window_size 500000000 ' \
                     f'--step_size 500000000 ' \
                     f'--resolution 25000 ' \
                     f'--batch_size {batch_size} ' \
                     f'--pad 2 --std 1 --rp 0.5 --tol 0.01 ' \
                     f'--min_cutoff {min_cutoff} ' \
                     f'--cpu_per_job {cpu_per_job}'
    print('Prepare imputation snakefiles for 25Kb resolution.')
    execute_command(impute_25k_cmd)

    # 100K
    this_out_dir = pathlib.Path(f'{impute_dir}/100K/')
    this_out_dir.mkdir(exist_ok=True)
    impute_100k_cmd = f'hicluster imputation ' \
                      f'--input_scool {scool_path_100k} ' \
                      f'--output_dir {this_out_dir} ' \
                      f'--chrom_size_path {chrom_size_path} ' \
                      f'--output_dist 500000000 ' \
                      f'--window_size 500000000 ' \
                      f'--step_size 500000000 ' \
                      f'--resolution 100000 ' \
                      f'--batch_size {batch_size} ' \
                      f'--pad 1 --std 1 --rp 0.5 --tol 0.01 ' \
                      f'--min_cutoff {min_cutoff} ' \
                      f'--cpu_per_job {cpu_per_job}'
    print('Prepare imputation snakefiles for 100Kb resolution.')
    execute_command(impute_100k_cmd)

    if scheduler == 'qsub':
        # prepare qsub
        qsub_dir = snakemake_dir / 'qsub'
        qsub_dir.mkdir(exist_ok=True)
        subprocess.run(f'cat {impute_dir}/*/snakemake_cmd.txt > {qsub_dir}/snakemake_cmd.txt',
                       shell=True, check=True)
        qsub_str = f'yap qsub --command_file_path {qsub_dir}/snakemake_cmd.txt ' \
                   f'--working_dir {qsub_dir} --project_name y{project_name} ' \
                   f'--total_cpu {total_cpu} --qsub_global_parms "-pe smp={cpu_per_job};-l h_vmem=5G"'
        with open(qsub_dir / 'qsub.sh', 'w') as f:
            f.write(qsub_str)
    elif scheduler == 'sbatch':
        # an important notice
        print(f'\n\n\nYou used sbatch as scheduler, remember to '
              f'1) provide the remote path of chrom_size_path; '
              f'2) transfer the {project_name}/scool directory to $SCRATCH/{project_name}/scool. '
              f'Otherwise the snakemake will raise error.')

        # prepare sbatch
        local_sbatch_dir = snakemake_dir / 'sbatch'
        remote_sbatch_dir = f'$SCRATCH/{project_name}/scool/snakemake/sbatch'
        local_sbatch_dir.mkdir(exist_ok=True)
        snake_command_files = impute_dir.glob('*/snakemake_cmd.txt')
        # concat all the snakemake commands and swap dir path
        with open(f'{local_sbatch_dir}/snakemake_cmd.txt', 'w') as out_f:
            for path in snake_command_files:
                with open(path) as f:
                    for line in f:
                        regex_txt = f'(?<= )\S+?(?=/{project_name}/scool/impute/)'
                        p = re.compile(regex_txt)
                        line = p.sub('$SCRATCH', line)
                        out_f.write(line)
        sbatch_str = f'yap sbatch --project_name {project_name} ' \
                     f'--command_file_path {remote_sbatch_dir}/snakemake_cmd.txt ' \
                     f'--working_dir {remote_sbatch_dir} --time_str {sbatch_time_str}'
        with open(local_sbatch_dir / 'sbatch.sh', 'w') as f:
            f.write(sbatch_str)

        # finally the path in snakefile need to be changed, go though each snakefile
        for path in impute_dir.glob('*/*/Snakefile'):
            modify_snakefile_for_sbatch(path, project_name)
    else:
        raise ValueError(f'scheduler can only be qsub or sbatch, got {scheduler}')
    return


def calculate_3c_datasets(output_dir, fasta_path, cpu=10):
    output_dir = pathlib.Path(output_dir).absolute()
    project_name = output_dir.name
    scool_dir = output_dir / 'scool'
    snakemake_dir = scool_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True, parents=True)
    impute_dir = scool_dir / 'impute'
    impute_dir.mkdir(exist_ok=True)
    dataset_dir = scool_dir / 'dataset'
    dataset_dir.mkdir(exist_ok=True)

    # validate imputation
    none_success = []
    for sub_path in impute_dir.iterdir():
        if sub_path.is_dir():
            for chunk_dir in sub_path.glob('chunk*'):
                flag = False
                for path in chunk_dir.iterdir():
                    if path.name == 'Success':
                        flag = True
                if not flag:
                    none_success.append(chunk_dir)
    if len(none_success) != 0:
        print('These imputation directories did not successfully executed, check and run imputation again:')
        for path in none_success:
            print(path)
        raise FileNotFoundError

    # Compartment calling command
    compartment_input_dir = impute_dir / '100K'
    compartment_cell_table = pd.Series({
        path.name.split('.')[0]: str(path)
        for path in compartment_input_dir.glob('*/*.cool')
    })
    compartment_cell_table_path = compartment_input_dir / 'cell_table.tsv'
    compartment_cell_table.to_csv(compartment_cell_table_path, sep='\t', header=None)
    cpg_path = compartment_input_dir / 'cpg_ratio.hdf'
    # compute CpG ratio first
    cpg_ratio_cmd = f'hicluster cpg-ratio --cell_url {compartment_cell_table.iloc[0]} ' \
                    f'--fasta_path {fasta_path} --hdf_output_path {cpg_path}'
    execute_command(cpg_ratio_cmd)
    compartment_cmd = f'hicluster compartment ' \
                      f'--cell_table_path {compartment_cell_table_path} ' \
                      f'--output_prefix {dataset_dir / project_name} ' \
                      f'--cpg_profile_path {cpg_path} ' \
                      f'--cpu {cpu}'

    # Domain calling command
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

    # Embedding command
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
    qsub_dir.mkdir(exist_ok=True)
    with open(qsub_dir / 'dataset_commands.txt', 'w') as f:
        f.write('\n'.join([compartment_cmd, domain_cmd, embedding_cmd]))
    qsub_str = f'yap qsub --command_file_path {qsub_dir}/snakemake_cmd.txt ' \
               f'--working_dir {qsub_dir} --project_name y{project_name} ' \
               f'--total_cpu {int(cpu * 3)} --qsub_global_parms "-pe smp={cpu};-l h_vmem=5G"'
    with open(qsub_dir / 'qsub.sh', 'w') as f:
        f.write(qsub_str)
    print('Run this command to generate datasets')
    print(f"sh {qsub_dir / 'qsub.sh'}")
    return
