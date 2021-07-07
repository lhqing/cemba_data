import pandas as pd
import pathlib
import yaml
from cooler.util import binnify, read_chromsizes
import numpy as np
import time
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
import subprocess
import pybedtools
from concurrent.futures import as_completed, ProcessPoolExecutor

import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])
DSS_TWO_GROUP_TEMPLATE = PACKAGE_DIR / 'dmr/dss/DSS.TwoGroup.SingleRegionDML.ipynb'


def prepare_snakemake(allc_table_path, output_dir, chrom_sizes_path, template_path, chroms=None, smoothing=True,
                      chunk_size=50000000):
    output_dir = pathlib.Path(output_dir).absolute()
    output_dir.mkdir(exist_ok=True)

    allc_table = pd.read_csv(allc_table_path, sep='\t')
    allc_table.columns = ['allc_path', 'sample', 'group']

    if allc_table['group'].unique().size != 2:
        raise ValueError(
            f"There must be two and only two different groups, got {allc_table['group'].unique().size}."
        )
    group1, group2 = allc_table['group'].unique()
    group1_allc = allc_table.loc[allc_table['group'] == group1,
                                 'allc_path'].tolist()
    group2_allc = allc_table.loc[allc_table['group'] == group2,
                                 'allc_path'].tolist()
    group1_id = allc_table.loc[allc_table['group'] == group1, 'sample'].tolist()
    group2_id = allc_table.loc[allc_table['group'] == group2, 'sample'].tolist()

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
        parameters = dict(region=region,
                          allc_paths=group1_allc + group2_allc,
                          group1=group1_id,
                          group2=group2_id,
                          smoothing=smoothing)
        with open(config_path, 'w') as f:
            f.write(yaml.dump(parameters))

    snakefile = f"""
regions = {regions}
rule main:
    input:
        expand('{{region}}.DSS.DML.hdf', region=regions)

rule papermill:
    input:
        nb='{template_path}',
        config='{{region}}.yaml'
    output:
        nb='{{region}}.ipynb',
        data='{{region}}.DSS.DML.hdf'
    shell:
        'papermill {{input.nb}} {{output.nb}} -f {{input.config}} && sleep 10'
"""
    snakefile_path = f'{output_dir}/Snakefile'
    with open(snakefile_path, 'w') as f:
        f.write(snakefile)
    return snakefile_path


def _parse_dml_ids(string):
    ids = np.array(string.split(',')).astype(int)
    start = min(ids)
    end = max(ids)
    result = pd.Series({
        'idx_start': start,
        'idx_end': end,
        'n_sig': ids.size,
        'total_dml': end - start + 1,
        'sig_ratio': ids.size / (end - start + 1)
    })
    return result


def call_dmr_single_chromosome(output_dir, chrom, p_threshold, min_cg, min_len,
                               sig_ratio, delta, chrom_sizes_path):
    # read DML for one chromosome
    print(f'Reading DML tables for {chrom}')
    dss_paths = list(pathlib.Path(output_dir).glob(f'{chrom}:*DSS.DML.hdf'))
    dmls = []
    for path in dss_paths:
        try:
            df = pd.read_hdf(path)
            if df.shape[0] > 0:
                dmls.append(df)
        except ValueError:
            # this happens when the HDF5 is empty
            pass
    dmls = pd.concat(dmls)
    dmls.sort_values('pos', inplace=True)
    dmls.dropna(subset=['pval'], inplace=True)
    dmls.reset_index(drop=True, inplace=True)

    print('Selecting significant DMLs to merge DMRs')
    # Select sig
    # recalculate FDR
    _, fdr, *_ = multipletests(dmls['pval'], method='fdr_bh')
    dmls['fdr'] = fdr
    # Calculate delta significance Pr(diff.abs() > delta)
    delta_p = 1 - norm.cdf(dmls['diff'].abs() - delta, scale=dmls['diff.se'])
    # select sig DML to merge DMR
    dmls['sig'] = (dmls['fdr'] < p_threshold) & (delta_p < p_threshold)
    dmls.to_hdf(f'{output_dir}/{chrom}.DML.hdf', key='data', format="table")
    sig_dmls = dmls[dmls['sig']]

    # Merge DMR
    print('Merge DMRs')
    dml_bed = sig_dmls[['chr', 'pos', 'pos']].reset_index()
    dml_bed.columns = ['id', 'chr', 'start', 'end']
    dml_bed['end'] += 1
    dml_bed = dml_bed.iloc[:, [1, 2, 3, 0]]
    dml_bed = pybedtools.BedTool.from_dataframe(dml_bed)
    try:
        dmr_bed = dml_bed.sort(g=chrom_sizes_path).merge(
            d=250, c='4', o='collapse').to_dataframe()
    except pybedtools.helpers.BEDToolsError:
        return False
    # Annotate DMR
    print('Annotating DMRs')
    name = dmr_bed.pop('name')
    dmr_bed = pd.concat([dmr_bed, name.apply(_parse_dml_ids)], axis=1)
    dmr_bed = dmr_bed.astype({
        'idx_start': int,
        'idx_end': int,
        'n_sig': int,
        'total_dml': int
    })

    def _region_stat(row):
        idx_start = row['idx_start']
        idx_end = row['idx_end']
        return dmls.iloc[idx_start:idx_end + 1].agg({
            'mu1': 'mean',
            'mu2': 'mean',
            'stat': 'sum'
        })

    dmr_stats = dmr_bed.apply(_region_stat, axis=1)
    dmr_bed = pd.concat([dmr_bed, dmr_stats], axis=1)
    dmr_bed['diff'] = dmr_bed['mu1'] - dmr_bed['mu2']
    dmr_bed['length'] = dmr_bed['end'] - dmr_bed['start']

    # final DMR filter
    judge = (dmr_bed['n_sig'] >= min_cg) & (
            dmr_bed['sig_ratio'] > sig_ratio) & (dmr_bed['length'] >= min_len) & (
                    dmr_bed['diff'].abs() > delta)
    dmr_bed['selected_dmr'] = judge
    dmr_bed.to_hdf(f'{output_dir}/{chrom}.DMR.hdf', key='data')
    return True


def run_dss_two_group(allc_table_path, output_dir, study_name, chrom_sizes_path, chroms=None, smoothing=True,
                      p_threshold=0.001, min_cg=1, min_len=1, sig_ratio=0.5, delta=0.1, cpu=10, save_dml=False,
                      template_path='default', chunk_size=50000000):
    # prepare template
    if template_path == 'default':
        template_path = DSS_TWO_GROUP_TEMPLATE

    # prepare snakemake
    output_dir = pathlib.Path(output_dir).absolute()
    this_study_dir = output_dir / f'{study_name}_DSS'
    prepare_snakemake(allc_table_path=allc_table_path,
                      output_dir=this_study_dir,
                      chrom_sizes_path=chrom_sizes_path,
                      chroms=chroms,
                      smoothing=smoothing,
                      template_path=template_path,
                      chunk_size=chunk_size)
    # the ipykernel package raise "zmq.error.ZMQError: Address already in use" due to parallel execution,
    # restart likely solve the problem.
    snakemake_cmd = f'snakemake -d {this_study_dir} --snakefile {this_study_dir}/Snakefile ' \
                    f'-j {cpu} --default-resources mem_mb=100 --resources mem_mb={int(5000 * cpu)} ' \
                    f'--restart-times 3'
    subprocess.run(snakemake_cmd, shell=True, check=True)

    with ProcessPoolExecutor(cpu) as exe:
        futures = {}
        for chrom in chroms:
            f = exe.submit(call_dmr_single_chromosome,
                           output_dir=this_study_dir,
                           chrom=chrom,
                           p_threshold=p_threshold,
                           min_cg=min_cg,
                           min_len=min_len,
                           sig_ratio=sig_ratio,
                           delta=delta,
                           chrom_sizes_path=chrom_sizes_path)
            futures[f] = chrom
            time.sleep(1)

        total_dmrs = []
        for future in as_completed(futures):
            chrom = futures[future]
            print(f'{chrom} finished')
            success_flag = future.result()
            if success_flag:
                total_dmrs.append(pd.read_hdf(f'{this_study_dir}/{chrom}.DMR.hdf'))
        if len(total_dmrs) != 0:
            total_dmrs = pd.concat(total_dmrs)
            total_dmrs.to_hdf(f'{this_study_dir}/{study_name}_DMR.hdf', key='data')

    for chrom in chroms:
        subprocess.run(f'rm -f {this_study_dir}/{chrom}:*.DSS.DML.hdf', shell=True)
        subprocess.run(f'rm -f {this_study_dir}/{chrom}:*.yaml', shell=True)
        subprocess.run(f'rm -f {this_study_dir}/{chrom}.DMR.hdf', shell=True)
        if not save_dml:
            subprocess.run(f'rm -f {this_study_dir}/{chrom}.DML.hdf', shell=True)
        subprocess.run(f'rm -rf {this_study_dir}/.snakemake {this_study_dir}/Snakefile', shell=True)
        subprocess.run(f'rm -f {this_study_dir}/*log {this_study_dir}/qsub.sh', shell=True)
    return
