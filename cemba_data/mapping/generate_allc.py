import pathlib
import warnings

import pandas as pd

from .utilities import parse_mc_pattern

from .utilities import get_mode
def generate_allc(
        output_dir,
        reference_fasta,
        num_upstr_bases=0,
        num_downstr_bases=2,
        compress_level=5):
    output_dir = pathlib.Path(output_dir).absolute()
    bam_batch = pd.read_csv(output_dir / 'snakemake/bismark_bam_list.txt',
                            header=None, index_col=0, squeeze=True)
    allc_dir = output_dir / 'allc'
    allc_dir.mkdir(exist_ok=True)
    mode = get_mode(output_dir)

    for batch_id, sub_series in bam_batch.groupby(bam_batch):
        bam_dict = {pathlib.Path(i).name.split('.')[0]: i
                    for i in sub_series.index}
        total_rules = ''
        allc_paths = []
        for i, (cell_id, bam_path) in enumerate(bam_dict.items()):
            uid = '-'.join(cell_id.split('-')[:-1])
            this_allc_dir = allc_dir / uid
            this_allc_dir.mkdir(exist_ok=True)
            allc_path = this_allc_dir / f'{cell_id}.allc.tsv.gz'
            stat_path = this_allc_dir / f'{cell_id}.allc.tsv.gz.count.csv'
            rule_template = f"""
rule allc_{i}:
    input:
        "{bam_path}"
    output:
        "{allc_path}"
    log:
        "{stat_path}"
    shell:
        'allcools bam-to-allc '
        '--bam_path {{input}} '
        '--reference_fasta {reference_fasta} '
        '--output_path {{output}} '
        '--cpu 1 '
        '--num_upstr_bases {num_upstr_bases} '
        '--num_downstr_bases {num_downstr_bases} '
        '--compress_level {compress_level} '
        '--save_count_df'

"""
            total_rules += rule_template
            allc_paths.append(str(allc_path))

        if mode == 'mc':
            include_str = f'include: "{output_dir}/snakemake/snakefile_bismark_mapping_{batch_id}"'
        elif mode == 'mct':
            include_str = f'include: "{output_dir}/snakemake/snakefile_select_dna_{batch_id}"'
        else:
            raise ValueError(f'Unknown mode {mode}')

        with open(output_dir / f'snakemake/snakefile_generate_allc_{batch_id}', 'w') as f:
            f.write(f"""
{include_str}

rule allc:
    input:
        {allc_paths}
    output:
        touch("{output_dir}/snakemake/generate_allc_done_{batch_id}")

{total_rules}
""")
    return


def _parse_allc_stat(stat_path, mc_patterns):
    try:
        report_df = pd.read_csv(stat_path, index_col=0)
    except pd.errors.EmptyDataError:
        report_df = pd.DataFrame([], columns=['mc', 'cov', 'mc_rate', 'genome_cov'])
    cell_id = pathlib.Path(stat_path).name.split('.')[0]
    report_df['cell_id'] = cell_id

    mc_rates = {}
    for pattern in mc_patterns:
        contexts = parse_mc_pattern(pattern)
        mc_sum = report_df.loc[report_df.index.isin(contexts), ['mc', 'cov']].sum()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mc_rates[f'{pattern}Rate'] = mc_sum['mc'] / mc_sum['cov']
    mc_rate = pd.Series(mc_rates, name=cell_id)
    return mc_rate


def generate_allc_stats(output_dir, patterns=('CHN', 'CGN', 'CCC')):
    output_dir = pathlib.Path(output_dir).absolute()
    allc_stats_dict = {p.name.split('.')[0]: p for p in output_dir.glob('allc/*/*count.csv')}
    # real all cell stats
    total_stats = []
    for cell_id, path in allc_stats_dict.items():
        allc_stat = pd.read_csv(path, index_col=0)
        allc_stat['cell_id'] = cell_id
        total_stats.append(allc_stat)
    total_stats = pd.concat(total_stats)

    # aggregate into patterns
    cell_records = []
    for pattern in patterns:
        contexts = parse_mc_pattern(pattern)
        pattern_stats = total_stats[total_stats.index.isin(contexts)]
        cell_sum = pattern_stats.groupby('cell_id')[['mc', 'cov']].sum()
        cell_rate = cell_sum['mc'] / cell_sum['cov']

        pattern_trans = {'CCC': 'mCCC',
                         'CGN': 'mCG',
                         'CHN': 'mCH'}
        _pattern = pattern_trans[pattern] if pattern in pattern_trans else pattern
        cell_rate.name = f'{_pattern}Rate'
        cell_records.append(cell_rate)
    final_df = pd.DataFrame(cell_records).T.reindex(allc_stats_dict.keys())
    final_df.index.name = 'cell_id'
    final_df.to_csv(output_dir / 'allc/allc_stats.csv')
    return
