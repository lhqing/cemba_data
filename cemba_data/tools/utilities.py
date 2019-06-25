import itertools
import functools
import collections
import pathlib
import numpy as np
import pandas as pd
from typing import Union, List

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


@functools.lru_cache(maxsize=100)
def parse_mc_pattern(pattern: str) -> set:
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


@functools.lru_cache(maxsize=10)
def parse_chrom_size(path, remove_chr_list=None):
    """
    support simple UCSC chrom size file, or .fai format (1st and 2nd columns same as chrom size file)

    return chrom:length dict
    """
    if remove_chr_list is None:
        remove_chr_list = []

    with open(path) as f:
        chrom_dict = collections.OrderedDict()
        for line in f:
            # *_ for other format like fadix file
            chrom, length, *_ = line.strip('\n').split('\t')
            if chrom in remove_chr_list:
                continue
            chrom_dict[chrom] = int(length)
    return chrom_dict


def genome_region_chunks(chrom_size_file: str,
                         bin_length: int=10000000,
                         combine_small: bool=True) -> List[str]:
    """
    Split the whole genome into bins, where each bin is {bin_length} bp. Used for tabix region query

    Parameters
    ----------
    chrom_size_file
        Path of UCSC genome size file
    bin_length
        length of each bin
    combine_small
        whether combine small regions into one record

    Returns
    -------
    list of records in tabix query format
    """
    chrom_size_dict = parse_chrom_size(chrom_size_file)

    cur_chrom_pos = 0
    records = []
    record_lengths = []
    for chrom, chrom_length in chrom_size_dict.items():
        while cur_chrom_pos + bin_length <= chrom_length:
            # tabix region is 1 based and inclusive
            records.append(f'{chrom}:{cur_chrom_pos}-{cur_chrom_pos+bin_length-1}')
            cur_chrom_pos += bin_length
            record_lengths.append(bin_length)
        else:
            records.append(f'{chrom}:{cur_chrom_pos}-{chrom_length}')
            cur_chrom_pos = 0
            record_lengths.append(chrom_length - cur_chrom_pos)

    # merge small records (when bin larger then chrom length)
    final_records = []
    if combine_small:
        temp_records = []
        cum_length = 0
        for record, record_length in zip(records, record_lengths):
            temp_records.append(record)
            cum_length += record_length
            if cum_length >= bin_length:
                final_records.append(' '.join(temp_records))
                temp_records = []
                cum_length = 0
        if len(temp_records) != 0:
            final_records.append(' '.join(temp_records))
    else:
        for record in records:
            final_records.append(record)
    return final_records


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


def parse_file_paths(input_file_paths: Union[str, list]) -> list:
    if isinstance(input_file_paths, list) and (len(input_file_paths) == 1):
        input_file_paths = input_file_paths[0]

    if isinstance(input_file_paths, str):
        if '*' in input_file_paths:
            import glob
            file_list = glob.glob(input_file_paths)
        else:
            file_list = []
            with open(input_file_paths) as f:
                for line in f:
                    file_list.append(line.strip('\n'))
        _file_list = file_list
    elif isinstance(input_file_paths, list):
        _file_list = input_file_paths
    else:
        raise TypeError('File paths input is neither str nor list.')

    final_file_list = []
    for path in _file_list:
        real_path = pathlib.Path(path).resolve()
        if not real_path.exists():
            raise FileNotFoundError(f'{path} provided do not exist.')
        final_file_list.append(str(real_path))
    return _file_list


def read_msg_bed(file_path):
    from pybedtools import BedTool

    df = pd.read_msgpack(file_path)
    bed = BedTool.from_dataframe(df)

    return bed
