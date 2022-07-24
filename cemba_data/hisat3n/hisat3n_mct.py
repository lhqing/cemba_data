import pathlib
from collections import defaultdict

import pandas as pd
import pysam


def _determine_reads_conversion(read):
    """
    Determine HISAT-3N read conversion type

    Parameters
    ----------
    read

    Returns
    -------
    return conversion string GA or CT
    """

    # For HISAT-3N, use YZ tag
    # YZ:A:<A>: The value + or – indicate the read is mapped to REF-3N (+) or REF-RC-3N (-).
    # actually, for snmC, the + means G to A conversion; the - means C to T conversion.
    # For HISAT2, just determine the read type
    # R1 is C to T conversion, R2 is G to A conversion
    try:
        yz_tag = read.get_tag('YZ')
        if yz_tag == '-':
            # G -> A
            return 'GA'
        elif yz_tag == '+':
            # C -> T
            return 'CT'
    except KeyError:
        if read.is_read1:
            return 'GA'
        elif read.is_read2:
            return 'CT'
        else:
            raise ValueError
    else:
        raise ValueError


def _complement(seq):
    d = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    try:
        return ''.join([d[b] for b in seq])
    except KeyError:
        return seq


def _determine_mch_context(conversion, ref_pos, ref_seq_dict, nome=False):
    """
    Determine whether the base is mCH context or not. If nome=True, skip GpC sites.

    Parameters
    ----------
    conversion :
        Type of conversion, GA or CT
    ref_pos :
        Position of the base in the reference
    ref_seq_dict :
        Dictionary of the reference sequence
    nome :
        If True, skip all the GpC context when calculating read mCH fraction,
        as it is subject to NOMe methylation

    Returns
    -------
    flag : bool
    """
    flag = False
    try:
        if conversion == 'GA':  # this means G to A conversion
            # the actual C is on the reverse strand
            ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos - 1]
            ref_context = _complement(ref_context)
        else:  # this means C to T conversion
            # the actual C is on the forward strand
            ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos + 1]

        if ref_context in {'CC', 'CT', 'CA'}:
            flag = True
            if nome:
                if conversion == 'GA':
                    if ref_seq_dict[ref_pos + 1] == 'G':
                        flag = False
                else:
                    if ref_seq_dict[ref_pos - 1] == 'G':
                        flag = False
    except KeyError:
        # the base is at boarder
        pass
    return flag


def _single_read_mch_level(read, nome=False, frac=False):
    """
    Determine the mCH fraction of single read

    Parameters
    ----------
    read
        Read object from bam file
    nome
        If True, skip all the GpC context as it is subject to methylation
    frac
        If True, return the mCH fraction of the read instead of mCH count

    Returns
    -------

    """
    read_conversion = _determine_reads_conversion(read)

    # ref seq is parsed based on read seq and MD tag, and do not depend on reverse or not
    ref_seq = read.get_reference_sequence().upper()
    ref_pos = read.get_reference_positions()
    # use dict instead of string is because ref_seq could contain blocks when skip happen
    ref_seq_dict = {pos: base for pos, base in zip(ref_pos, ref_seq)}
    read_seq = read.seq.upper()

    # only count mCH
    mch = 0
    cov = 0
    other_snp = 0
    # read.get_aligned_pairs always return small to large position
    for read_pos, ref_pos, ref_base in read.get_aligned_pairs(
            matches_only=True, with_seq=True):
        read_base = read_seq[read_pos]
        ref_read_pair = ref_base + read_base
        ref_read_pair = ref_read_pair.upper()
        is_mch_context = _determine_mch_context(read_conversion, ref_pos, ref_seq_dict, nome=nome)
        if is_mch_context:
            if ref_read_pair in ['CC', 'GG']:  # means unconverted and methylated
                cov += 1
                mch += 1
            elif ref_read_pair in ['CT', 'GA']:  # means converted and un-methylated
                cov += 1
            else:
                # other kinds of SNPs, do not count to cov
                other_snp += 1
                pass
    if frac:
        read_mch_frac = (mch / cov) if cov > 0 else 0
        return read_mch_frac, cov, other_snp
    else:
        return mch, cov, other_snp


def select_mct_reads(input_bam,
                     output_bam,
                     mode,
                     mc_rate_max_threshold=None,
                     mc_rate_min_threshold=None,
                     cov_min_threshold=3,
                     nome=False):
    """
    Select DNA reads with mCH fraction <= mc_rate_max_threshold and
    coverage >= cov_min_threshold from mCT bam file

    Parameters
    ----------
    input_bam :
        Input bam file
    output_bam :
        Output bam file
    mode :
        Mode of reads selection, 'dna' or 'rna'
    mc_rate_max_threshold :
        Maximum mCH fraction cutoff, applied only when mode is 'dna'
    mc_rate_min_threshold :
        Minimum mCH fraction cutoff, applied only when mode is 'rna'
    cov_min_threshold :
        Minimum coverage cutoff
    nome :
        If True, skip all the GpC context when calculating read mCH fraction,
        as it is subject to NOMe methylation

    Returns
    -------

    """
    if mode == 'dna':
        if mc_rate_max_threshold is None:
            raise ValueError('mc_rate_max_threshold is required when mode is dna')
    elif mode == 'rna':
        if mc_rate_min_threshold is None:
            raise ValueError('mc_rate_min_threshold is required when mode is rna')
    else:
        raise ValueError('mode should be dna or rna')

    read_profile_dict = defaultdict(int)
    # init dict to make sure the series has something
    read_profile_dict[(50, 50)] = 0

    with pysam.AlignmentFile(input_bam) as f:
        with pysam.AlignmentFile(output_bam, header=f.header,
                                 mode='wb') as out_f:
            for read in f:
                mc_frac, cov, _ = _single_read_mch_level(read, nome=nome, frac=True)
                read_profile_dict[(int(100 * mc_frac), cov)] += 1

                # split reads
                if mode == 'dna':
                    # determine dna reads
                    if (mc_frac > mc_rate_max_threshold) or (cov < cov_min_threshold):
                        continue
                else:
                    # determine rna reads
                    if (mc_frac < mc_rate_min_threshold) or (cov < cov_min_threshold):
                        continue
                out_f.write(read)

    with open(str(output_bam)[:-4] + '.reads_mch_frac.csv', 'w') as stat_f:
        # save parameters, so the stat parser can read from it
        stat_f.write(f'#mode={mode}\n')
        stat_f.write(f'#mc_rate_max_threshold={mc_rate_max_threshold}\n')
        stat_f.write(f'#mc_rate_min_threshold={mc_rate_min_threshold}\n')
        stat_f.write(f'#cov_min_threshold={cov_min_threshold}\n')

        stat_f.write('mc_frac,cov,count\n')
        for (mc_frac, cov), count in read_profile_dict.items():
            stat_f.write(f'{mc_frac},{cov},{count}\n')
    return


def aggregate_feature_counts():
    cell_datas = []

    save_info = True
    for path in pathlib.Path('rna_bam/').glob('*.feature_count.tsv'):
        table = pd.read_csv(path, comment='#', sep='\t', index_col=0)
        cell_data = table.iloc[:, -1].squeeze()
        cell_data.name = cell_data.name.split(':')[-1]
        cell_datas.append(cell_data)
        if save_info:
            table.iloc[:, :-1].to_csv('featureCounts.gene_info.csv.gz')
            save_info = False

    pd.DataFrame(cell_datas).to_csv('featureCounts.data.csv.gz')
    return
