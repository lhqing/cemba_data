import pathlib
import subprocess
from collections import defaultdict

import pandas as pd
import pysam

from .utilities import get_configuration

REVERSE_READ_MCH_CONTEXT = {'CA', 'CC', 'CT'}
FORWARD_READ_MCH_CONTEXT = {'AG', 'TG', 'GG'}


def reverse_complement(seq):
    _seq = ''
    base_map = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C'}
    for char in seq.upper()[::-1]:
        _seq += base_map[char]
    return _seq


def single_read_mch_level(read):
    # ref seq is parsed based on read seq and MD tag, and do not depend on reverse or not
    ref_seq = read.get_reference_sequence().upper()
    # use dict instead of string is because ref_seq could contain blocks when skip happen
    ref_seq_dict = {pos: base for pos, base in zip(read.get_reference_positions(), ref_seq)}
    read_seq = read.seq.upper()

    # only count mCH
    mch = 0
    cov = 0
    other_snp = 0
    if read.is_reverse:  # read in reverse strand
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos + 1]
                if ref_context not in REVERSE_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue
            if ref_read_pair == 'CC':  # C to C means unconverted and methylated
                cov += 1
                mch += 1
            elif ref_read_pair == 'CT':  # C to T means converted and un-methylated
                cov += 1
            else:
                # other kinds of SNPs, do not count to cov
                other_snp += 1
                pass
    else:  # read in forward strand
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos - 1] + ref_seq_dict[ref_pos]
                if ref_context not in FORWARD_READ_MCH_CONTEXT:
                    continue
            except KeyError:
                # ref_seq_dict KeyError means position is on border or not continuous, skip that
                continue
            if ref_read_pair == 'GG':  # G to G means unconverted and methylated
                cov += 1
                mch += 1
            elif ref_read_pair == 'GA':  # G to A means converted and un-methylated
                cov += 1
            else:
                # other kinds of SNPs, do not count to cov
                other_snp += 1
                pass
    read_mch_rate = (mch / cov) if cov > 0 else 0
    return read_mch_rate, cov, other_snp


def select_rna_reads(input_bam,
                     output_bam,
                     mc_rate_min_threshold=0.9,
                     cov_min_threshold=5,
                     remove_input=True):
    read_profile_dict = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as bam:
        with pysam.AlignmentFile(output_bam, header=bam.header, mode='w') as out_bam:
            for read in bam:
                read_mch_rate, cov, other_snp = single_read_mch_level(read)
                read_profile_dict[(int(100 * read_mch_rate), cov)] += 1

                # split reads
                if (read_mch_rate < mc_rate_min_threshold) or (cov < cov_min_threshold):
                    continue
                out_bam.write(read)

    read_profile = pd.Series(read_profile_dict)
    read_profile.index.name = ['mc_rate', 'cov']
    read_profile.to_csv(str(output_bam) + '.reads_profile.csv', header=True)
    if remove_input:
        subprocess.run(['rm', '-f', input_bam])
    return


def prepare_select_rna_reads(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    if isinstance(config, str):
        config = get_configuration(config)

    bismark_records = pd.read_csv(output_dir / 'star_bam_qc.records.csv',
                                  index_col=['uid', 'index_name'],
                                  squeeze=True)
    mc_rate_min_threshold = config['RNAReadsFilter']['mc_rate_min_threshold']
    cov_min_threshold = config['RNAReadsFilter']['cov_min_threshold']
    remove_input = config['RNAReadsFilter']['remove_input']

    # process bam
    records = []
    command_list = []
    for (uid, index_name), bismark_bam_path in bismark_records.iteritems():
        # file path
        output_bam = bismark_bam_path[:-3] + 'rna_reads.bam'
        # command
        # TODO change this to PE mapping, right now only map R1

        keep_input_str = '--remove_input' if remove_input else ''
        command = f'yap-internal select-rna-reads ' \
                  f'--input_bam {bismark_bam_path} ' \
                  f'--output_bam {output_bam} ' \
                  f'--mc_rate_min_threshold {mc_rate_min_threshold} ' \
                  f'--cov_min_threshold {cov_min_threshold} ' \
                  f'{keep_input_str}'
        records.append([uid, index_name, output_bam])
        command_list.append(command)

    with open(output_dir / 'select_rna_reads.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'bam_path'])
    record_df.to_csv(output_dir / 'select_rna_reads.records.csv', index=None)
    return record_df, command_list


def summarize_select_rna_reads(output_dir, config):
    bam_dir = pathlib.Path(output_dir)
    output_path = bam_dir / 'select_rna_reads.stats.csv'
    if output_path.exists():
        return str(output_path)
    config = get_configuration(config)

    records = []
    select_dna_reads_stat_list = list(bam_dir.glob('*.reads_profile.csv'))
    for path in select_dna_reads_stat_list:
        try:
            _df = pd.read_csv(path)
            if _df.shape[0] == 0:
                subprocess.run(['rm', '-f', path])
                continue
            _df.columns = ['mch_rate', 'ch_cov', 'read_count']
            report_series = pd.Series({col_name: '|'.join(col.astype(str))
                                       for col_name, col in _df.iteritems()})
        except pd.errors.EmptyDataError:
            # means the bam file is empty
            subprocess.run(['rm', '-f', path])
            continue

        *uid, index_name, suffix = path.name.split('-')
        uid = '-'.join(uid)
        read_type = suffix.split('.')[0]
        report_series['uid'] = uid
        report_series['index_name'] = index_name
        report_series['read_type'] = read_type
        report_series['mc_rate_min_threshold'] = config['RNAReadsFilter']['mc_rate_min_threshold']
        report_series['cov_min_threshold'] = config['RNAReadsFilter']['cov_min_threshold']
        records.append(report_series)
        subprocess.run(['rm', '-f', path])
    total_stats_df = pd.DataFrame(records)
    total_stats_df.to_csv(output_path, index=None)
    return str(output_path)
