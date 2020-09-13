import pathlib

import pandas as pd

from cemba_data.utilities import get_barcode_version

from cemba_data.utilities import get_configuration

def add_v1_plateinfo(total_stats):
    total_plate_pos_records = {}
    for cell_id in total_stats.index:
        *uid, index_name = cell_id.split('-')
        uid = '-'.join(uid)
        record_dict = {}
        plate1, plate2, pos96 = uid.split('-')
        record_dict['plate1'] = plate1
        record_dict['plate2'] = plate2
        record_dict['Pos96'] = pos96
        index_name = index_name.lower()
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
        total_plate_pos_records[cell_id] = record_dict

    plate_info_df = pd.DataFrame(total_plate_pos_records).T
    total_data = pd.concat([total_stats, plate_info_df], axis=1)
    return total_data


def add_v2_plateinfo(total_stats):
    total_plate_pos_records = {}
    for cell_id in total_stats.index:
        *uid, index_name = cell_id.split('-')
        uid = '-'.join(uid)
        record_dict = {}
        plate, multiplex_group, primer_name = uid.split('-')
        record_dict['multiplex_group'] = int(multiplex_group)
        # primer name is a Pos384 of the illumina primer,
        # but this pos has nothing to do with cell Pos384.
        # Cell Pos384 determine by random primer pos (index_name)
        record_dict['primer_name'] = primer_name
        # V2 doesn't cross multiplex plate, real_plate is plate
        record_dict['real_plate'] = plate

        # 384 pos
        record_dict['Pos384'] = index_name
        record_dict['Col384'] = int(record_dict['Pos384'][1:]) - 1
        record_dict['Row384'] = ord(record_dict['Pos384'][0]) - 65  # convert A-H to 0-8
        total_plate_pos_records[cell_id] = record_dict
    plate_info_df = pd.DataFrame(total_plate_pos_records).T
    total_data = pd.concat([total_stats, plate_info_df], axis=1)
    return total_data


def basic_summary(output_dir, patterns=('CHN', 'CGN', 'CCC')):
    barcode_version = get_barcode_version(output_dir)

    output_dir = pathlib.Path(output_dir).absolute()
    demultiplex_stat = pd.read_csv(output_dir / 'fastq/demultiplex.stats.csv',
                                   index_col='cell_id')

    from .bismark_mapping import bismark_mapping_stats
    bismark_mapping_stats(output_dir)
    mapping_stat = pd.read_csv(output_dir / 'bam/bismark_mapping_stats.csv',
                               index_col='cell_id')

    from .generate_allc import generate_allc_stats
    generate_allc_stats(output_dir, patterns=patterns)
    allc_stat = pd.read_csv(output_dir / 'allc/allc_stats.csv',
                            index_col='cell_id')

    total_stats = pd.concat([demultiplex_stat, mapping_stat, allc_stat],
                            axis=1)

    # add plate info
    if barcode_version == 'V1':
        total_stats = add_v1_plateinfo(total_stats)
    elif barcode_version == 'V2':
        total_stats = add_v2_plateinfo(total_stats)
    else:
        raise ValueError
    return total_stats


def snmc_summary(output_dir, patterns=('CHN', 'CGN', 'CCC')):
    total_stats = basic_summary(output_dir, patterns=patterns)
    total_stats.to_csv(f'{output_dir}/MappingSummary.csv.gz')


def snmct_summary(output_dir, patterns=('CHN', 'CGN', 'CCC')):
    basic_stats = basic_summary(output_dir, patterns=patterns)
    # summarize rna and dna reads selection
    from .star_mapping import summary_rna_mapping
    summary_rna_mapping(output_dir)
    rna_stats = pd.read_csv(f'{output_dir}/rna_bam/star_mapping_stats.csv', index_col=0)

    from .mct_bismark_bam_filter import summarize_select_dna_reads
    summarize_select_dna_reads(output_dir,
                               mc_rate_max_threshold=0.5,
                               cov_min_threshold=3)
    dna_stats = pd.read_csv(f'{output_dir}/dna_bam/select_dna_reads.stats.csv', index_col=0)

    final_stats = pd.concat([basic_stats, dna_stats, rna_stats], axis=1)
    final_stats.to_csv(f'{output_dir}/MappingSummary.csv.gz')


def summary(output_dir):
    config = get_configuration(config_path=f'{output_dir}/snakemake/mapping_config.ini')
    mode = config['mode'].lower()

    if mode == 'mc':
        snmc_summary(output_dir, patterns=config['mc_stat_feature'])
    elif mode == 'mct':
        snmct_summary(output_dir)

    # TODO add cleaning function
    return
