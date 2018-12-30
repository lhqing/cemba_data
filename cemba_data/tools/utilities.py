import itertools
import collections
import numpy as np

IUPAC_TABLE = {
    'A': 'A',
    'T': 'T',
    'C': 'C',
    'G': 'G',
    'R': 'AG',
    'Y': 'CT',
    'S': 'GC',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ATC',
    'V': 'ACG',
    'N': 'ATCG'
}


def parse_mc_pattern(pattern):
    """
    parse mC context pattern
    """
    # IUPAC DNA abbr. table
    all_pos_list = []
    pattern = pattern.upper()
    for base in pattern:
        try:
            all_pos_list.append(IUPAC_TABLE[base])
        except KeyError:
            raise KeyError(f'Base {base} is not in IUPAC table.')
    context_set = set([''.join(i) for i in itertools.product(*all_pos_list)])
    return context_set


def parse_chrom_size(path, remove_chr_list=None):
    """
    return chrom:length dict
    """
    if remove_chr_list is None:
        remove_chr_list = []
    with open(path) as f:
        chrom_dict = collections.OrderedDict()
        for line in f:
            chrom, length = line.strip('\n').split('\t')
            if chrom in remove_chr_list:
                continue
            chrom_dict[chrom] = int(length)
    return chrom_dict


def get_mean_var(X):
    # - using sklearn.StandardScaler throws an error related to
    #   int to long trafo for very large matrices
    # - using X.multiply is slower
    mean = X.mean(axis=0)
    # scanpy deal with both sparse and full matrix, here only support full
    # if issparse(X):
    #     mean_sq = X.multiply(X).mean(axis=0)
    #     mean = mean.A1
    #     mean_sq = mean_sq.A1
    # else:
    #     mean_sq = np.multiply(X, X).mean(axis=0)
    mean_sq = np.multiply(X, X).mean(axis=0)
    # enforce R convention (unbiased estimator) for variance
    var = (mean_sq - mean ** 2) * (X.shape[0] / (X.shape[0] - 1))
    return mean, var
