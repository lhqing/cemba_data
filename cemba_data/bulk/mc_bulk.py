import pandas as pd
import pathlib
import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def prepare_mc_bulk(allc_table,
                    output_dir,
                    bigwig_contexts,
                    chrom_size_path,
                    extract_mcg=False,
                    bigwig_bin_size=50,
                    merge_allc_cpu=8,
                    mcg_context='CGN'):
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
    extract_mcg
        Whether run the step to extract mCG sites from the merged ALLC.
        If your input ALLC only contains mCG sites, this can be skipped.
        Otherwise, this needs to be done before running the CG-DMR calling.
    bigwig_contexts
        mC contexts for generating the bigwig tracks
    chrom_size_path
        Path of the chromosome size file path
    bigwig_bin_size
        Bin size used to generate bigwig
    merge_allc_cpu
        Number of CPU to use in individual merge-allc job
    mcg_context
        mC context for extract_mcg step, only relevant when extract_mcg=True

    Returns
    -------

    """
    snakemake_template_path = PACKAGE_DIR / 'bulk/Snakefile_template/mc_bulk.Snakefile'
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

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
    allc_dir = output_dir / 'allc'
    allc_dir.mkdir(exist_ok=True)

    for group, paths in allc_dict.items():
        allc_list_path = allc_dir / f'{group}.allc_paths.txt'
        with open(allc_list_path, 'w') as f:
            f.write('\n'.join(paths))

    snakemake_parameters = f"""
merge_allc_cpu = {merge_allc_cpu}
extract_mcg = {extract_mcg}
mcg_context = '{mcg_context}'
bigwig_contexts = {bigwig_contexts}
bigwig_bin_size = {bigwig_bin_size}
chrom_size_path = '{chrom_size_path}'
groups = {list(allc_dict.keys())}

"""

    with open(snakemake_template_path) as f:
        snakemake_template = f.read()
    snakemake_str = snakemake_parameters + snakemake_template

    with open(output_dir / 'Snakefile', 'w') as f:
        f.write(snakemake_str)
    return
