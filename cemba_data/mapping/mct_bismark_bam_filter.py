import subprocess
from collections import defaultdict

import pandas as pd
from ALLCools._open import open_bam

from .utilities import get_bam_header_str

METHYLATED_CHAR = 'H'
UNMETHYLATED_CHAR = 'h'


def read_mc_level(bismark_tag):
    m_c = sum([bismark_tag.count(c) for c in METHYLATED_CHAR])
    normal_c = sum([bismark_tag.count(c) for c in UNMETHYLATED_CHAR])
    total_c = m_c + normal_c
    if total_c == 0:
        return 0, 0
    read_mc_rate = int(100 * (m_c / total_c))
    return read_mc_rate, total_c


def filter_bismark_reads_mc_level(input_bam,
                                  output_bam,
                                  mc_rate_max_threshold=0.5,
                                  cov_min_threshold=5,
                                  remove_input=True):
    bam_header = get_bam_header_str(input_bam)
    read_profile_dict = defaultdict(int)
    with open_bam(input_bam, include_header=False) as f, open_bam(output_bam, 'w') as out_f:
        out_f.write(bam_header)
        for line in f:
            bismark_tag = line.split('\t')[-3]
            mc_rate, cov = read_mc_level(bismark_tag)
            read_profile_dict[(mc_rate, cov)] += 1

            # split reads
            if (mc_rate > mc_rate_max_threshold) or (cov < cov_min_threshold):
                continue
            out_f.write(line)
    read_profile = pd.Series(read_profile_dict)
    read_profile.index.name = ['mc_rate', 'cov']
    read_profile.to_csv(str(output_bam) + '.reads_profile.csv', header=True)
    if remove_input:
        subprocess.run(['rm', '-f', input_bam])
    return
