from collections import defaultdict

import pysam

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
                     cov_min_threshold=3):
    read_profile_dict = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as bam:
        with pysam.AlignmentFile(output_bam, header=bam.header, mode='wb') as out_bam:
            for read in bam:
                read_mch_rate, cov, other_snp = single_read_mch_level(read)
                read_profile_dict[(int(100 * read_mch_rate), cov)] += 1

                # split reads
                if (read_mch_rate < mc_rate_min_threshold) or (cov < cov_min_threshold):
                    continue
                out_bam.write(read)

    with open(str(output_bam) + '.reads_profile.csv', 'w') as stat_f:
        stat_f.write('mc_rate,cov,count\n')
        for (mc_rate, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_rate},{cov},{count}\n')
    return
