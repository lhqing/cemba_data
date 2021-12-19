from collections import defaultdict
import re
import pandas as pd
import pysam


def determine_reads_conversion(read):
    """
    Determine HISAT-3N read conversion level

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
    yz_tag = read.get_tag('YZ')
    if yz_tag == '-':
        # G -> A
        return 'GA'
    elif yz_tag == '+':
        # C -> T
        return 'CT'
    else:
        raise ValueError


def complement(seq):
    d = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    try:
        return ''.join([d[b] for b in seq])
    except KeyError:
        return seq


def determine_mch_context(conversion, ref_pos, ref_seq_dict, nome=False):
    """
    Determine whether the base is mCH context or not. If nome=True, skip GpC sites.
    """
    flag = False
    try:
        if conversion == 'GA':  # this means G to A conversion
            # the actual C is on the reverse strand
            ref_context = ref_seq_dict[ref_pos] + ref_seq_dict[ref_pos - 1]
            ref_context = complement(ref_context)
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


def single_read_mch_level(read, nome=False, frac=False):
    """
    Determine the mCH fraction of single read

    Parameters
    ----------
    read
    nome
        If True, skip all the GpC context as it is subject to methylation
    frac

    Returns
    -------

    """
    read_conversion = determine_reads_conversion(read)

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
        is_mch_context = determine_mch_context(read_conversion, ref_pos, ref_seq_dict, nome=nome)
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
        return mch, cov, other_snp
    else:
        read_mch_frac = (mch / cov) if cov > 0 else 0
        return read_mch_frac, cov, other_snp
