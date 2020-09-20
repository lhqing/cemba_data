from collections import defaultdict

import pysam

METHYLATED_CHAR = 'H'
UNMETHYLATED_CHAR = 'h'


def read_mc_level(bismark_tag):
    m_c = bismark_tag.count(METHYLATED_CHAR)
    normal_c = bismark_tag.count(UNMETHYLATED_CHAR)
    total_c = m_c + normal_c
    if total_c == 0:
        return 0, 0
    else:
        read_mc_rate = m_c / total_c
        return read_mc_rate, total_c


def select_dna_reads(input_bam,
                     output_bam,
                     mc_rate_max_threshold=0.5,
                     cov_min_threshold=3):
    read_profile_dict = defaultdict(int)
    # init dict to make sure the series has something
    read_profile_dict[(50, 50)] = 0
    with pysam.AlignmentFile(input_bam) as f:
        with pysam.AlignmentFile(output_bam, header=f.header, mode='wb') as out_f:
            for read in f:
                bismark_tag = read.get_tag('XM')
                mc_rate, cov = read_mc_level(bismark_tag)
                read_profile_dict[(int(100 * mc_rate), cov)] += 1

                # split reads
                if (mc_rate > mc_rate_max_threshold) or (cov < cov_min_threshold):
                    continue
                out_f.write(read)
    with open(str(output_bam) + '.reads_profile.csv', 'w') as stat_f:
        stat_f.write('mc_rate,cov,count\n')
        for (mc_rate, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_rate},{cov},{count}\n')
    return
