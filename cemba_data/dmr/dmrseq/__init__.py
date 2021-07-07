import pandas as pd
import pathlib
import yaml
from cooler.util import binnify, read_chromsizes
import subprocess
from statsmodels.stats.multitest import multipletests

import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])
DMRSEQ_TEMPLATE = PACKAGE_DIR / 'dmr/dmrseq/DMRseq.ipynb'


def prepare_snakemake(allc_table_path,
                      output_dir,
                      chrom_sizes_path,
                      template_path,
                      chroms=None,
                      test_covariate='group',
                      match_covariate=None,
                      adjust_covariate=None,
                      cutoff=0.1,
                      min_num_region=3,
                      smooth=True,
                      bp_span=1000,
                      min_in_span=30,
                      max_gap_smooth=2500,
                      max_gap=1000,
                      verbose=True,
                      max_perms=10,
                      stat="stat",
                      block=False,
                      block_size=5000,
                      chrs_per_chunk=1,
                      cpu=40,
                      chunk_size=5000000000):
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

    chrom_sizes = read_chromsizes(chrom_sizes_path).reindex(chroms)
    if chroms is None:
        chroms = chrom_sizes.index.tolist()
    bins = binnify(chrom_sizes.loc[chroms], binsize=chunk_size)
    regions = []
    for _, (chrom, start, end) in bins.iterrows():
        region = f'{chrom}:{start}-{end}'
        regions.append(region)

    for region in regions:
        config_path = f'{output_dir}/{region}.yaml'
        parameters = {
            'region': region,
            'allc_table_path': allc_table_path,
            'test_covariate': test_covariate,
            'match_covariate': match_covariate,
            'adjust_covariate': adjust_covariate,
            'cutoff': cutoff,
            'min_num_region': min_num_region,
            'smooth': smooth,
            'bp_span': bp_span,
            'min_in_span': min_in_span,
            'max_gap_smooth': max_gap_smooth,
            'max_gap': max_gap,
            'verbose': verbose,
            'max_perms': max_perms,
            'stat': stat,
            'block': block,
            'block_size': block_size,
            'chrs_per_chunk': chrs_per_chunk,
            'cpu': cpu
        }
        with open(config_path, 'w') as f:
            f.write(yaml.dump(parameters))

    snakefile = f"""
regions = {regions}
rule main:
    input:
        expand('{{region}}.DMR.hdf', region=regions)

rule papermill:
    input:
        nb='{template_path}',
        config='{{region}}.yaml'
    output:
        nb='{{region}}.ipynb',
        data='{{region}}.DMR.hdf'
    threads:
        1 #{cpu}
    shell:
        'papermill {{input.nb}} {{output.nb}} -f {{input.config}} && sleep 10'
"""
    snakefile_path = f'{output_dir}/Snakefile'
    with open(snakefile_path, 'w') as f:
        f.write(snakefile)
    return regions


def run_dmrseq(allc_table_path,
               output_dir,
               study_name,
               chrom_sizes_path,
               chroms=None,
               test_covariate='group',
               match_covariate=None,
               adjust_covariate=None,
               cutoff=0.1,
               min_num_region=3,
               smooth=True,
               bp_span=1000,
               min_in_span=30,
               max_gap_smooth=2500,
               max_gap=1000,
               verbose=True,
               max_perms=10,
               stat="stat",
               block=False,
               block_size=5000,
               chrs_per_chunk=1,
               cpu=40,
               template_path='default',
               chunk_size=50000000):
    # prepare template
    if template_path == 'default':
        template_path = DMRSEQ_TEMPLATE
    allc_table_path = str(pathlib.Path(allc_table_path).absolute())

    # prepare snakemake
    output_dir = pathlib.Path(output_dir).absolute()
    this_study_dir = output_dir / f'{study_name}_DMRseq'
    regions = prepare_snakemake(allc_table_path=allc_table_path,
                                output_dir=this_study_dir,
                                chrom_sizes_path=chrom_sizes_path,
                                template_path=template_path,
                                chroms=chroms,
                                test_covariate=test_covariate,
                                match_covariate=match_covariate,
                                adjust_covariate=adjust_covariate,
                                cutoff=cutoff,
                                min_num_region=min_num_region,
                                smooth=smooth,
                                bp_span=bp_span,
                                min_in_span=min_in_span,
                                max_gap_smooth=max_gap_smooth,
                                max_gap=max_gap,
                                verbose=verbose,
                                max_perms=max_perms,
                                stat=stat,
                                block=block,
                                block_size=block_size,
                                chrs_per_chunk=chrs_per_chunk,
                                cpu=cpu,
                                chunk_size=chunk_size)
    # the ipykernel package raise "zmq.error.ZMQError: Address already in use" due to parallel execution,
    # restart likely solve the problem.
    snakemake_cmd = f'snakemake -d {this_study_dir} --snakefile {this_study_dir}/Snakefile ' \
                    f'-j {cpu} --default-resources mem_mb=100 --resources mem_mb={int(5000 * cpu)} ' \
                    f'--restart-times 3'
    subprocess.run(snakemake_cmd, shell=True, check=True)

    total_dmrs = []
    for region in regions:
        try:
            total_dmrs.append(pd.read_hdf(f'{this_study_dir}/{region}.DMR.hdf'))
        except ValueError:
            pass
    # happens when there is no sites to test in dmrseq, and output file is empty
    total_dmrs = pd.concat(total_dmrs)
    # recalculate FDR
    # some pval in the dmrseq result are nan
    total_dmrs = total_dmrs.loc[~total_dmrs['pval'].isna()].copy()
    _, fdr, *_ = multipletests(total_dmrs['pval'], method='fdr_bh')
    total_dmrs['qval'] = fdr
    total_dmrs.to_hdf(f'{this_study_dir}/{study_name}_DMR.hdf', key='data', format="table")

    for region in regions:
        subprocess.run(f'rm -f {this_study_dir}/{region}.DMR.hdf',
                       shell=True)
        subprocess.run(f'rm -f {this_study_dir}/*.yaml', shell=True)
        subprocess.run(
            f'rm -rf {this_study_dir}/.snakemake {this_study_dir}/Snakefile',
            shell=True)
    return
