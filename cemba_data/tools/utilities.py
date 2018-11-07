import itertools
import collections


def parse_mc_pattern(pattern):
    """
    parse mC context pattern
    :param pattern: length and C position should match ALLC
    :return: context set
    """
    all_pos_list = []
    for base in pattern:
        if base in 'hH':
            all_pos_list.append('ATC')
        elif base in '*-_Nn':
            all_pos_list.append('ATCG')
        else:
            all_pos_list.append(base.upper())
    context_set = set([''.join(i) for i in itertools.product(*all_pos_list)])
    return context_set


def parse_chrom_size(path, remove_chr_list=None):
    """
    return chrom:length dict
    :param path: ucsc chrom.size file
    :return: chrom:length dict
    """
    if remove_chr_list is None:
        remove_chr_list = []
    with open(path) as f:
        chrom_dict = collections.OrderedDict()
        for line in f:
            chrom, length = line.strip('\n').split('\t')
            if add_chr:
                if ('chr' != chrom[:3]) and add_chr:
                    chrom = 'chr' + chrom
            if chrom in remove_chr_list:
                continue
            chrom_dict[chrom] = int(length)
    return chrom_dict
