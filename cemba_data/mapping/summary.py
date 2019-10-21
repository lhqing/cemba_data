import pathlib

import numpy as np
import pandas as pd

from cemba_data.mapping import get_configuration


def transform_demultiplex_stats(stats_df_dict):
    try:
        stats_df = stats_df_dict['demultiplex']
    except KeyError:
        return None

    # demultiplex stat, from cutadapt demultiplex step
    demultiplex_result = stats_df.groupby(['uid', 'index_name']) \
        .sum()[['TotalPair', 'Trimmed']]
    demultiplex_result.rename(columns={'TotalPair': 'MultiplexReadsTotal',
                                       'Trimmed': 'IndexReadsTotal'},
                              inplace=True)
    demultiplex_result['MultiplexReadsTotal'] *= 2  # R1 and R2
    demultiplex_result['IndexReadsTotal'] *= 2  # R1 and R2
    demultiplex_result['IndexReadsRatio'] = demultiplex_result['IndexReadsTotal'] / demultiplex_result[
        'MultiplexReadsTotal'] * 100
    return demultiplex_result


def transform_fastq_qc_stats(stats_df_dict, stats_final_dict):
    stats_df = stats_df_dict['fastq_qc']
    demultiplex_result = stats_final_dict['demultiplex']
    # fastq trim stat
    fastq_trim_result = stats_df.groupby(['uid', 'index_name']).sum()[
        ['in_bp', 'out_bp', 'out_reads', 'qualtrim_bp', 'too_short', 'w/adapters']]
    fastq_trim_result.rename(columns={'in_bp': 'IndexBpTotal',
                                      'out_bp': 'IndexTrimedBpTotal',
                                      'out_reads': 'IndexTrimedReadsTotal',
                                      'qualtrim_bp': 'ReadsQualTrimBpTotal',
                                      'too_short': 'ReadsLengthFilterTotal',
                                      'w/adapters': 'ReadsWithAdapterTotal'}, inplace=True)
    fastq_trim_result['TrimedReadsAveLength'] = fastq_trim_result['IndexTrimedBpTotal'] / fastq_trim_result[
        'IndexTrimedReadsTotal']
    if demultiplex_result is not None:
        fastq_trim_result['IndexTrimedReadsRatio'] = fastq_trim_result['IndexTrimedReadsTotal'] / demultiplex_result[
            'IndexReadsTotal']
    else:
        fastq_trim_result['IndexTrimedReadsRatio'] = np.NaN
    return fastq_trim_result


def transform_bismark_mapping_stats(stats_df_dict):
    stats_df = stats_df_dict['bismark_mapping']
    # bismark_mapping stat
    bismark_r1 = stats_df[stats_df['read_type'].apply(lambda i: '1' in i)] \
        .set_index(['uid', 'index_name'])[['CTOB', 'CTOT', 'mapping_rate', 'total_c',
                                           'total_reads', 'unique_map', 'unmap', 'ununique_map']]
    bismark_r1.rename(columns={'mapping_rate': 'R1MappedRatio',  # this mapping rate is unique mapping rate
                               'total_c': 'R1TotalC',
                               'total_reads': 'R1TrimmedReads',
                               'unique_map': 'R1UniqueMappedReads',
                               'unmap': 'R1UnmappedReads',
                               'ununique_map': 'R1UnuniqueMappedReads'},
                      inplace=True)
    bismark_r2 = stats_df[stats_df['read_type'].apply(lambda i: '2' in i)] \
        .set_index(['uid', 'index_name'])[['OB', 'OT', 'mapping_rate', 'total_c',
                                           'total_reads', 'unique_map', 'unmap', 'ununique_map']]
    bismark_r2.rename(columns={'mapping_rate': 'R2MappedRatio',
                               'total_c': 'R2TotalC',
                               'total_reads': 'R2TrimmedReads',
                               'unique_map': 'R2UniqueMappedReads',
                               'unmap': 'R2UnmappedReads',
                               'ununique_map': 'R2UnuniqueMappedReads'},
                      inplace=True)
    bismark_result = pd.concat([bismark_r1, bismark_r2], sort=True, axis=1)
    bismark_result['TotalUniqueMappedReads'] = \
        bismark_result['R1UniqueMappedReads'] + bismark_result['R2UniqueMappedReads']
    bismark_result['TotalMappedRatio'] = bismark_result['TotalUniqueMappedReads'] / bismark_result[[
        'R1TrimmedReads', 'R2TrimmedReads']].sum(axis=1)
    return bismark_result


def transform_bismark_bam_qc_stats(stats_df_dict):
    stats_df = stats_df_dict['bismark_bam_qc']
    bam_result = stats_df.groupby(['uid', 'index_name']) \
        .sum()[['UNPAIRED_READ_DUPLICATES', 'UNPAIRED_READS_EXAMINED']] \
        .rename(columns={'UNPAIRED_READ_DUPLICATES': 'DupReads',
                         'UNPAIRED_READS_EXAMINED': 'TotalReadsBeforeDeDup'})
    bam_result['DeduppedReads'] = bam_result['TotalReadsBeforeDeDup'] - bam_result['DupReads']
    bam_result['DeduppedRatio'] = bam_result['DeduppedReads'] / bam_result['TotalReadsBeforeDeDup']
    return bam_result


def _stats_df_profile(stats_df):
    cov_df = stats_df \
        .set_index(['uid', 'index_name', 'mc_context'])['cov'] \
        .unstack('mc_context') \
        .fillna(0).astype(int)
    ccc_cov = cov_df[f'CCC']
    cov_df = cov_df.groupby(lambda i: i[:2], axis=1) \
        .sum()[['CA', 'CC', 'CG', 'CT']] \
        .rename(columns={c: c + '_Cov'
                         for c in ['CA', 'CC', 'CG', 'CT']})

    mc_df = stats_df \
        .set_index(['uid', 'index_name', 'mc_context'])['mc'] \
        .unstack('mc_context') \
        .fillna(0).astype(int)
    ccc_mc = mc_df['CCC']
    mc_df = mc_df.groupby(lambda i: i[:2], axis=1) \
        .sum()[['CA', 'CC', 'CG', 'CT']] \
        .rename(columns={c: c + '_Mc'
                         for c in ['CA', 'CC', 'CG', 'CT']})

    # mCCC
    mc_df['CCC_Mc'] = ccc_mc
    mc_df['CCC_Cov'] = ccc_cov
    mc_df['CCC_Rate'] = ccc_mc / ccc_cov

    # mCG
    mc_df['CG_Rate'] = mc_df['CG_Mc'] / cov_df['CG_Cov']
    mc_df['CG_RateAdj'] = (mc_df['CG_Rate'] - mc_df['CCC_Rate']) / (1 - mc_df['CCC_Rate'])

    # mCH
    cov_df[f'CH_Cov'] = cov_df[['CA_Cov', 'CC_Cov', 'CT_Cov']].sum(axis=1)
    mc_df['CH_Mc'] = mc_df[['CA_Mc', 'CC_Mc', 'CT_Mc']].sum(axis=1)
    mc_df['CH_Rate'] = mc_df['CH_Mc'] / cov_df['CH_Cov']
    mc_df['CH_RateAdj'] = (mc_df['CH_Rate'] - mc_df['CCC_Rate']) / (1 - mc_df['CCC_Rate'])

    # mCY
    cov_df[f'CY_Cov'] = cov_df[['CC_Cov', 'CT_Cov']].sum(axis=1)
    mc_df['CY_Mc'] = mc_df[['CC_Mc', 'CT_Mc']].sum(axis=1)
    mc_df['CY_Rate'] = mc_df['CY_Mc'] / cov_df['CY_Cov']
    mc_df['CY_RateAdj'] = (mc_df['CY_Rate'] - mc_df['CCC_Rate']) / (1 - mc_df['CCC_Rate'])
    return mc_df


def transform_allc_stats(stats_df_dict, nome=False):
    stats_df = stats_df_dict['generate_allc']
    stats_df.set_index(stats_df.columns[0], inplcae=True)
    stats_df.index.name = 'mc_context'
    stats_df.reset_index(inplace=True)
    genome_cov = stats_df.groupby(['uid', 'index_name'])['genome_cov'].apply(lambda i: i.iloc[0])

    if nome:
        nome_stats_df = stats_df[stats_df.mc_context.str.startswith('G')].copy()  # select GCN records
        nome_stats_df['mc_context'] = nome_stats_df['mc_context'].str[1:]  # remove the first base
        nome_results = _stats_df_profile(nome_stats_df)

        stats_df = stats_df[~stats_df.mc_context.str.startswith('G')].copy()  # select HCN records
        stats_df['mc_context'] = stats_df['mc_context'].str[1:]  # remove the first base
        stats_df = stats_df.groupby(['uid', 'index_name', 'mc_context']).sum()[['mc', 'cov']].reset_index()
        mc_results = _stats_df_profile(stats_df)

        # select useful col
        mc_results = mc_results[['CCC_Rate',
                                 'CG_Rate', 'CG_RateAdj',
                                 'CH_Rate', 'CH_RateAdj',
                                 'CY_Rate', 'CY_RateAdj']].copy()
        mc_results.columns = mc_results.columns.map(lambda i: f'H{i}')
        nome_results = nome_results[['CG_Rate',
                                     'CH_Rate',
                                     'CY_Rate']].copy()
        nome_results.columns = nome_results.columns.map(lambda i: f'G{i}')
        mc_results = pd.concat([nome_results, mc_results], axis=1)
        mc_results['GCH_RateAdj'] = (mc_results['GCH_Rate'] - mc_results['HCCC_Rate']) / (1 - mc_results['HCCC_Rate'])
        mc_results['GCY_RateAdj'] = (mc_results['GCY_Rate'] - mc_results['HCCC_Rate']) / (1 - mc_results['HCCC_Rate'])

    else:
        mc_results = _stats_df_profile(stats_df)
        mc_results = mc_results[['CCC_Rate',
                                 'CG_Rate', 'CG_RateAdj',
                                 'CH_Rate', 'CH_RateAdj']].copy()

    mc_results['genome_cov'] = genome_cov
    return mc_results


def transform_select_reads_stats(stats_df_dict, read_kind):
    stats_df = stats_df_dict[f'select_{read_kind.lower()}_reads'].set_index(['uid', 'index_name'])

    read_type_dfs = []
    for read_type, sub_df in stats_df.groupby('read_type'):
        total_records = {}
        for (uid, index_name), row in sub_df.iterrows():
            row_dict = {}
            read_type = row['read_type']

            row_df = pd.DataFrame([count_row
                                   for count_row in zip(
                    *row[['mch_rate', 'ch_cov', 'read_count']].str.split('|'))],
                                  columns=['mch_rate', 'ch_cov', 'read_count']).astype(int)
            total_reads = row_df['read_count'].sum()
            if read_kind.lower() == 'dna':
                pass_reads = row_df[(row_df['mch_rate'].astype(float) < row['mc_rate_max_threshold']) & (
                        row_df['ch_cov'].astype(float) > row['cov_min_threshold'])]['read_count'].sum()
            elif read_kind.lower() == 'rna':
                pass_reads = row_df[(row_df['mch_rate'].astype(float) > row['mc_rate_min_threshold']) & (
                        row_df['ch_cov'].astype(float) > row['cov_min_threshold'])]['read_count'].sum()
            else:
                raise ValueError(f'Unknown read_kind {read_kind}')
            row_dict['bismark_bam_final_reads'] = total_reads
            row_dict['selected_dna_reads'] = pass_reads
            total_records[(uid, index_name)] = row_dict
        reads_count_df = pd.DataFrame(total_records).T
        sub_df = pd.concat([sub_df, reads_count_df], sort=True, axis=1)
        sub_df.columns = sub_df.columns.map(lambda i: f'{read_kind}_{read_type}_{i}')
        read_type_dfs.append(sub_df)
    total_results = pd.concat(read_type_dfs, axis=1, sort=True)
    return total_results


def transform_star_mapping_stats(stats_df_dict):
    stats_df = stats_df_dict['star_mapping'].set_index(['uid', 'index_name'])
    read_type_dfs = []
    for read_type, sub_df in stats_df.groupby('read_type'):
        use_df = sub_df[['Uniquely mapped reads number', 'Number of splices: Total']]
        use_df.columns = use_df.columns.map(lambda i: f'{read_type}_{i}')
        read_type_dfs.append(use_df)
    total_results = pd.concat(read_type_dfs, axis=1)
    return total_results


POSSIBLE_STATS_NAME = [
    'demultiplex',
    'fastq_qc',
    'bismark_mapping',
    'bismark_bam_qc',
    'select_dna_reads',
    'generate_allc',
    'star_mapping',
    'select_rna_reads'
]


def aggregate_all_summary(output_dir, mct=False, nome=False):
    stats_df_dict = {}
    for path in pathlib.Path(output_dir).glob('**/*stats.csv'):
        name = path.name.split('.')[0]
        if name not in POSSIBLE_STATS_NAME:
            print(f'Found unknown name of stats file: {name}')
            continue
        stats_df_dict[name] = pd.read_csv(path)
    stats_final_dict = dict()

    stats_final_dict['demultiplex'] = transform_demultiplex_stats(stats_df_dict)
    stats_final_dict['fastq_qc'] = transform_fastq_qc_stats(stats_df_dict, stats_final_dict)
    stats_final_dict['bismark_mapping'] = transform_bismark_mapping_stats(stats_df_dict)
    stats_final_dict['bismark_bam_qc'] = transform_bismark_bam_qc_stats(stats_df_dict)
    stats_final_dict['generate_allc'] = transform_allc_stats(stats_df_dict, nome=nome)

    if mct:
        stats_final_dict['select_dna_reads'] = transform_select_reads_stats(stats_df_dict, 'dna')
        stats_final_dict['select_rna_reads'] = transform_select_reads_stats(stats_df_dict, 'rna')
        stats_final_dict['star_mapping'] = transform_star_mapping_stats(stats_df_dict)

    # concat total meta
    total_meta = pd.concat(stats_final_dict.values(),
                           sort=True, axis=1).dropna()
    return total_meta


def mapping_summary(output_dir):
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bismark_bam_dir = output_dir / 'bismark_bam'
    allc_dir = output_dir / 'allc'

    config_path = list(output_dir.glob('*ini'))[0]
    config = get_configuration(config_path)
    mct = 'mct' in config['mode']['mode'].lower()
    nome = 'nome' in config['mode']['mode'].lower()

    # summarize fastq
    from .demultiplex import summarize_demultiplex
    summarize_demultiplex(fastq_dir, config=config_path)
    from .fastq_qc import summarize_fastq_qc
    summarize_fastq_qc(fastq_dir)

    # summarize mapping
    from .bismark_mapping import summarize_bismark_mapping
    summarize_bismark_mapping(bismark_bam_dir)
    from .bismark_bam_qc import summarize_bismark_bam_qc
    summarize_bismark_bam_qc(bismark_bam_dir)

    # mCT specific
    if mct:
        from .mct_bismark_bam_filter import summarize_select_dna_reads
        summarize_select_dna_reads(bismark_bam_dir, config)

        from .star_mapping import summarize_star_mapping
        star_bam_dir = output_dir / 'star_bam'
        summarize_star_mapping(star_bam_dir)
        from .mct_star_bam_filter import summarize_select_rna_reads
        summarize_select_rna_reads(star_bam_dir, config)

    # summarize allc
    from .generate_allc import summarize_generate_allc
    summarize_generate_allc(allc_dir)

    # aggregate all files
    total_summary = aggregate_all_summary(output_dir, mct=mct, nome=nome)

    summary_path = output_dir / 'MappingSummary.csv.gz'
    total_summary.to_csv(summary_path)

    random_index_version = config['multiplexIndex']['barcode_version'].upper()
    if random_index_version == 'V1':
        total_summary_with_plate_info = add_v1_plateinfo(total_summary)
    elif random_index_version == 'V2':
        total_summary_with_plate_info = add_v2_plateinfo(total_summary)
    else:
        raise ValueError(f'Unknown version name {random_index_version} in multiplexIndex section of the config file.')

    total_summary_with_plate_info.to_csv(summary_path)
    return summary_path


def add_v1_plateinfo(total_summary):
    total_plate_pos_records = {}
    for _, (uid, index_name) in total_summary.reset_index()[['uid', 'index_name']].iterrows():
        record_dict = {}
        plate1, plate2, pos96 = uid.split('-')
        record_dict['plate1'] = plate1
        record_dict['plate2'] = plate2
        record_dict['Pos96'] = pos96
        index_name = index_name.lower()
        record_dict['index_name'] = index_name
        # judge real plate
        if index_name in ['ad001', 'ad002', 'ad004', 'ad006']:
            record_dict['real_plate'] = plate1
        else:
            record_dict['real_plate'] = plate2
        # 96 pos
        record_dict['Col96'] = int(record_dict['Pos96'][1:]) - 1
        record_dict['Row96'] = ord(record_dict['Pos96'][0]) - 65  # convert A-H to 0-8
        # 384 pos
        ad_index_384_dict = {
            'ad001': (0, 0), 'ad002': (0, 1),
            'ad004': (1, 0), 'ad006': (1, 1),
            'ad007': (0, 0), 'ad008': (0, 1),
            'ad010': (1, 0), 'ad012': (1, 1)}
        record_dict['Col384'] = 2 * record_dict['Col96'] + ad_index_384_dict[index_name][0]
        record_dict['Row384'] = 2 * record_dict['Row96'] + ad_index_384_dict[index_name][1]
        record_dict['Pos384'] = chr(record_dict['Row384'] + 65) + str(record_dict['Col384'])
        total_plate_pos_records[(uid, index_name)] = record_dict
    plate_info_df = pd.DataFrame(total_plate_pos_records).T
    plate_info_df.index = plate_info_df.index.set_names(['uid', 'index_name'])

    # make index_name lower to agree with plate_info_df
    total_summary.reset_index(inplace=True)
    total_summary['index_name'] = total_summary['index_name'].str.lower()
    total_summary = total_summary.set_index(['uid', 'index_name'])

    total_data = pd.concat([total_summary, plate_info_df], axis=1)
    return total_data


def add_v2_plateinfo(total_summary):
    total_plate_pos_records = {}
    for _, (uid, index_name) in total_summary.reset_index()[['uid', 'index_name']].iterrows():
        record_dict = {}
        plate, multiplex_group, primer_name = uid.split('-')
        record_dict['multiplex_group'] = int(multiplex_group)
        # primer name is a Pos384 of the illumina primer,
        # but this pos has nothing to do with cell Pos384.
        # Cell Pos384 determine by random primer pos (index_name)
        record_dict['primer_name'] = primer_name
        record_dict['index_name'] = index_name
        # V2 doesn't cross multiplex plate, real_plate is plate
        record_dict['real_plate'] = plate

        # 384 pos
        record_dict['Pos384'] = index_name
        record_dict['Col384'] = int(record_dict['Pos384'][1:]) - 1
        record_dict['Row384'] = ord(record_dict['Pos384'][0]) - 65  # convert A-H to 0-8
        total_plate_pos_records[(uid, index_name)] = record_dict
    plate_info_df = pd.DataFrame(total_plate_pos_records).T
    total_data = pd.concat([total_summary, plate_info_df], axis=1)
    return total_data
