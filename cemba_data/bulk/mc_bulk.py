import pandas as pd
import pathlib
import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def prepare_mc_bulk(allc_table,
                    output_dir,
                    chrom_size_path,
                    mch_context='CHN',
                    mcg_context='CGN',
                    bigwig_mch_bin_size=50,
                    bigwig_mcg_bin_size=1,
                    cpu_per_job=12,
                    total_cpu=60):
    """
    Prepare the snakefile for merging single-cell ALLC files into pseudo-bulk

    Parameters
    ----------
    allc_table
        Path of the allc table. The allc table is a two column tsv file.
        The first columns is the absolute ALLC file paths;
        the second column is the group name of each file.
    output_dir
        Path of the output directory, will be created if not exist.
    mch_context
        mCH contexts for generating the bigwig tracks
    mcg_context
        mCG contexts for generating the bigwig tracks and merge strand
    chrom_size_path
        Path of the chromosome size file path
    bigwig_mch_bin_size
        Bin size used to generate mCH bigwig
    bigwig_mcg_bin_size
        Bin size used to generate mCG bigwig
    cpu_per_job
        Number of CPUs to use in individual merge-allc job
    total_cpu
        Number of CPUs to use in total

    Returns
    -------

    """
    snakemake_template_path = PACKAGE_DIR / 'bulk/Snakefile_template/mc_bulk.Snakefile'
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

    merge_allc_cpu = int(cpu_per_job / 1.1)
    total_mem_mb = cpu_per_job * 5000

    # prepare ALLC path dict
    # allc_path to group
    allc_path = pd.read_csv(allc_table, sep='\t', index_col=0, header=None, squeeze=True)
    file_not_exist = allc_path[allc_path.index.map(lambda i: not pathlib.Path(i).exists())]
    if file_not_exist.size != 0:
        path_str = "\n".join(file_not_exist.index.tolist())
        raise FileNotFoundError(f'{file_not_exist.size} files do not exist:'
                                f'\n{path_str}')
    allc_dict = {group: paths.index.tolist() for group, paths in allc_path.groupby(allc_path)}

    # Prepare Snakefile
    snakemake_cmds = []
    for group, paths in allc_dict.items():
        # each group has a separate snakemake file
        group_dir = output_dir / group
        group_dir.mkdir(exist_ok=True)
        allc_list_path = group_dir / f'{group}.allc_paths.txt'
        with open(allc_list_path, 'w') as f:
            f.write('\n'.join(paths))
        snakemake_parameters = f"""
merge_allc_cpu = {merge_allc_cpu}
mch_context = '{mch_context}'
mcg_context = '{mcg_context}'
bigwig_mch_bin_size = {bigwig_mch_bin_size}
bigwig_mch_bin_size = {bigwig_mcg_bin_size}
chrom_size_path = '{chrom_size_path}'
group = '{group}'

"""
        with open(snakemake_template_path) as f:
            snakemake_template = f.read()
        snakemake_str = snakemake_parameters + snakemake_template
        with open(group_dir / f'Snakefile', 'w') as f:
            f.write(snakemake_str)
        snakemake_cmd = f'snakemake ' \
                        f'-d {group_dir.absolute()} ' \
                        f'--snakefile {group_dir.absolute()}/Snakefile ' \
                        f'-j {cpu_per_job} ' \
                        f'--default-resources mem_mb=100 ' \
                        f'--resources mem_mb={total_mem_mb} ' \
                        f'--rerun-incomplete'
        snakemake_cmds.append(snakemake_cmd)

    qsub_dir = output_dir / 'qsub'
    qsub_dir.mkdir(exist_ok=True)
    with open(qsub_dir / 'snakemake_cmds.txt', 'w') as f:
        f.write('\n'.join(snakemake_cmds))
    with open(qsub_dir / 'qsub.sh', 'w') as f:
        qsub_str = f"""
yap qsub \
--command_file_path {qsub_dir / 'snakemake_cmds.txt'} \
--working_dir {qsub_dir} \
--project_name merge \
--total_cpu {total_cpu} \
--qsub_global_parms "-pe smp={cpu_per_job};-l h_vmem=5G"
"""
        f.write(qsub_str)
        print(f'Execute this command to start pipeline:\nnohup sh {qsub_dir / "qsub.sh"} &')
    return
