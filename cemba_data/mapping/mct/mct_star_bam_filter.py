from collections import defaultdict
import re
import pandas as pd
import pysam

REVERSE_READ_MCH_CONTEXT = {'CA', 'CC', 'CT'}
FORWARD_READ_MCH_CONTEXT = {'AG', 'TG', 'GG'}


def single_read_mch_level(read, nome=False, frac=False):
    # ref seq is parsed based on read seq and MD tag, and do not depend on reverse or not
    ref_seq = read.get_reference_sequence().upper()
    # use dict instead of string is because ref_seq could contain blocks when skip happen
    ref_seq_dict = {
        pos: base
        for pos, base in zip(read.get_reference_positions(), ref_seq)
    }
    read_seq = read.seq.upper()

    # only count mCH
    mch = 0
    cov = 0
    other_snp = 0
    if read.is_reverse:  # read in reverse strand
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos + 1]
                if nome:
                    if ref_seq_dict[ref_pos - 1] == 'G':
                        continue
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
        for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
                matches_only=True, with_seq=True):
            read_base = read_seq[read_pos]
            ref_read_pair = ref_base + read_base
            try:
                ref_context = ref_seq_dict[ref_pos - 1] + ref_seq_dict[ref_pos]
                if nome:
                    if ref_seq_dict[ref_pos + 1] == 'C':
                        continue
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
    if frac:
        return mch, cov, other_snp
    else:
        read_mch_frac = (mch / cov) if cov > 0 else 0
        return read_mch_frac, cov, other_snp


def select_rna_reads_normal(input_bam,
                            output_bam,
                            mc_rate_min_threshold=0.9,
                            cov_min_threshold=3,
                            nome=False):
    read_profile_dict = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as bam:
        with pysam.AlignmentFile(output_bam, header=bam.header,
                                 mode='wb') as out_bam:
            for read in bam:
                read_mch_rate, cov, other_snp = single_read_mch_level(read, nome=nome)
                read_profile_dict[(int(100 * read_mch_rate), cov)] += 1

                # split reads
                if (read_mch_rate <
                    mc_rate_min_threshold) or (cov < cov_min_threshold):
                    continue
                out_bam.write(read)

    with open(str(output_bam) + '.reads_profile.csv', 'w') as stat_f:
        stat_f.write('mc_frac,cov,count\n')
        for (mc_rate, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_rate},{cov},{count}\n')
    return


def select_rna_reads_split_reads(input_bam,
                                 output_bam,
                                 mc_rate_min_threshold=0.9,
                                 cov_min_threshold=3,
                                 nome=False):
    splited_read_name_pattern = re.compile('.+-[lrm]$')

    # first pass: determine read methylation level
    read_level_mcs = defaultdict(int)
    read_level_covs = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as f:
        for read in f:
            mc, cov, other_snp = single_read_mch_level(read, frac=True, nome=nome)
            read_name = read.qname
            if splited_read_name_pattern.search(read_name):
                read_level_mcs[read_name[:-2]] += mc
                read_level_covs[read_name[:-2]] += cov
            else:
                read_level_mcs[read_name] += mc
                read_level_covs[read_name] += cov
    read_level_data = pd.DataFrame({
        'mc': read_level_mcs,
        'cov': read_level_covs
    })
    read_level_data['mc_frac'] = read_level_data['mc'] / (read_level_data['cov'] +
                                                       0.001)
    read_level_data['mc_frac'] = (read_level_data['mc_frac'] * 100).astype(int)
    profile = read_level_data.groupby('mc_frac')['cov'].value_counts()
    profile.name = 'count'
    profile = profile.reset_index()
    profile.to_csv(f'{output_bam}.reads_profile.csv', index=None)

    # filter reads
    use_reads = read_level_data[
        (read_level_data['mc_frac'] > mc_rate_min_threshold)
        & (read_level_data['cov'] >= cov_min_threshold)].index.tolist()
    use_reads = set(use_reads)
    del read_level_data

    # second pass: write passed reads
    with pysam.AlignmentFile(input_bam) as f:
        with pysam.AlignmentFile(output_bam, header=f.header,
                                 mode='wb') as out_f:
            for read in f:
                read_name = read.qname
                if (read_name in use_reads) or (read_name[:-2] in use_reads):
                    # read name or read name without suffix
                    out_f.write(read)
    return


def select_rna_reads(input_bam,
                     output_bam,
                     mc_rate_min_threshold=0.5,
                     cov_min_threshold=3,
                     nome=False,
                     assay_type='mc'):
    if assay_type == 'mc':
        select_rna_reads_normal(input_bam,
                                output_bam,
                                mc_rate_min_threshold=mc_rate_min_threshold,
                                cov_min_threshold=cov_min_threshold,
                                nome=nome)
    elif assay_type == 'm3c':
        select_rna_reads_split_reads(input_bam,
                                     output_bam,
                                     mc_rate_min_threshold=mc_rate_min_threshold,
                                     cov_min_threshold=cov_min_threshold,
                                     nome=nome)
    else:
        raise ValueError(f'Unknown assay_type {assay_type}.')
    return
