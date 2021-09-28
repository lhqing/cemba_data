from collections import defaultdict
import re
import pysam
import pandas as pd


def read_mc_level(read, frac=True, nome=False):
    bismark_tag = read.get_tag('XM')
    if nome:
        m_c = 0
        normal_c = 0
        seq = read.seq.upper()
        read_length = len(seq)
        for pos, xm_base in enumerate(bismark_tag):
            if xm_base in '.ZzUu':
                # skip unrelated base (.), CpG (Zz), CpUnknown (Uu)
                continue
            # Skip GpC
            try:
                if read.is_reverse:
                    if (pos == read_length) or (read.seq[pos + 1] == 'C'):
                        continue
                else:
                    if (pos == 0) or (read.seq[pos - 1] == 'G'):
                        continue
            except IndexError:
                # start or end of the read
                continue
            if xm_base in 'xh':
                normal_c += 1
            elif xm_base in 'XH':
                m_c += 1
            else:
                pass
    else:
        m_c = bismark_tag.count('X') + bismark_tag.count('H')
        normal_c = bismark_tag.count('x') + bismark_tag.count('h')

    total_c = m_c + normal_c
    if total_c == 0:
        return 0, 0
    else:
        if frac:
            read_mc_rate = m_c / total_c
            return read_mc_rate, total_c
        else:
            return m_c, total_c


def select_dna_reads_normal(input_bam,
                            output_bam,
                            mc_rate_max_threshold=0.5,
                            cov_min_threshold=3,
                            nome=False):
    read_profile_dict = defaultdict(int)
    # init dict to make sure the series has something
    read_profile_dict[(50, 50)] = 0
    with pysam.AlignmentFile(input_bam) as f:
        with pysam.AlignmentFile(output_bam, header=f.header,
                                 mode='wb') as out_f:
            for read in f:
                mc_frac, cov = read_mc_level(read, nome=nome)
                read_profile_dict[(int(100 * mc_frac), cov)] += 1

                # split reads
                if (mc_frac > mc_rate_max_threshold) or (cov <
                                                         cov_min_threshold):
                    continue
                out_f.write(read)
    with open(str(output_bam) + '.reads_profile.csv', 'w') as stat_f:
        stat_f.write('mc_frac,cov,count\n')
        for (mc_frac, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_frac},{cov},{count}\n')
    return


def select_dna_reads_split_reads(input_bam,
                                 output_bam,
                                 mc_rate_max_threshold=0.5,
                                 cov_min_threshold=3,
                                 nome=False):
    splited_read_name_pattern = re.compile('.+-[lrm]$')

    # first pass: determine read methylation level
    read_level_mcs = defaultdict(int)
    read_level_covs = defaultdict(int)
    with pysam.AlignmentFile(input_bam) as f:
        for read in f:
            mc, cov = read_mc_level(read, frac=False, nome=nome)
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
    if read_level_data.shape[0] == 0:
        # in case there is no read at all:
        with open(f'{output_bam}.reads_profile.csv', 'w') as f:
            f.write('mc_frac,cov,count\n')
            f.write('0,1,0\n')
    else:
        profile = read_level_data.groupby('mc_frac')['cov'].value_counts()
        profile.name = 'count'
        profile = profile.reset_index()
        profile.to_csv(f'{output_bam}.reads_profile.csv', index=None)

    # filter reads
    use_reads = read_level_data[
        (read_level_data['mc_frac'] < mc_rate_max_threshold)
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


def select_dna_reads(input_bam,
                     output_bam,
                     mc_rate_max_threshold=0.5,
                     cov_min_threshold=3,
                     nome=False,
                     assay_type='mc'):
    if assay_type == 'mc':
        select_dna_reads_normal(input_bam,
                                output_bam,
                                mc_rate_max_threshold=mc_rate_max_threshold,
                                cov_min_threshold=cov_min_threshold,
                                nome=nome)
    elif assay_type == 'm3c':
        select_dna_reads_split_reads(input_bam,
                                     output_bam,
                                     mc_rate_max_threshold=mc_rate_max_threshold,
                                     cov_min_threshold=cov_min_threshold,
                                     nome=nome)
    else:
        raise ValueError(f'Unknown assay_type {assay_type}.')
    return
