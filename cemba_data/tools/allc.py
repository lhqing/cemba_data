from .utilities import *
import gzip
from functools import partial
from pybedtools import BedTool, cleanup
from subprocess import run
from collections import defaultdict
import pandas as pd
import numpy as np
import os
from .methylpy_utilities import merge_allc_files


def _split_to_chrom_bed(allc_path, context_pattern, genome_size_path,
                        out_path_prefix, max_cov_cutoff=None):
    """
    Split ALLC into bed format, chrom column contain "chr".
    :param allc_path: Single ALLC file path
    :param context_pattern: comma separate context patterns or list
    :param out_path_prefix: Single output prefix
    :param max_cov_cutoff: 2 for single cell, None for bulk or merged allc
    :param genome_size_path: UCSC chrom size file path
    :return: Path dict for out put files
    """
    chrom_set = set(parse_chrom_size(genome_size_path).keys())

    # deal with some old allc file that don't have chr in chrom name
    ref_chrom_have_chr = False
    for chrom in chrom_set:
        if chrom.startswith('chr'):
            ref_chrom_have_chr = True
        # if any chrom in the genome_size_path have chr, treat the whole reference as ref_chrom_have_chr
    # whether add chr or not depending on this, judged later in the first line
    need_to_add_chr = False

    # prepare context
    if isinstance(context_pattern, str):
        context_pattern = context_pattern.split(',')
    pattern_dict = {c: parse_mc_pattern(c) for c in context_pattern}

    # prepare out path
    path_dict = {(c, chrom): out_path_prefix + f'.{c}.{chrom}.bed'
                 for c in context_pattern
                 for chrom in chrom_set}

    # open all paths
    open_func = partial(open, mode='a')
    # open func for allc:
    if '.gz' in allc_path[-3:]:
        allc_open_func = partial(gzip.open, mode='rt')
    else:
        allc_open_func = partial(open, mode='r')

    # split ALLC
    first = True
    cur_chrom = None
    handle_dict = None
    with allc_open_func(allc_path) as allc:
        for line in allc:
            if first:
                chrom = line.split('\t')[0]
                # judge if the first line have chr or not, if not,
                # but ref_chrom_have_chr is true, add chr for every line
                if ref_chrom_have_chr and not chrom.startswith('chr'):
                    need_to_add_chr = True
                    chrom = 'chr' + chrom
                if chrom not in chrom_set:
                    continue
                first = False
                cur_chrom = chrom
                handle_dict = {c: open_func(path_dict[(c, cur_chrom)]) for c in context_pattern}
            ll = line.split('\t')
            # filter max cov (for single cell data)
            if (max_cov_cutoff is not None) and (int(ll[5]) > max_cov_cutoff):
                continue
            # judge chrom
            chrom = ll[0]
            if need_to_add_chr:
                chrom = 'chr' + chrom
            if chrom not in chrom_set:
                continue
            if chrom != cur_chrom:
                cur_chrom = chrom
                for handle in handle_dict.values():
                    handle.close()
                handle_dict = {c: open_func(path_dict[(c, cur_chrom)]) for c in context_pattern}
            # bed format [chrom, start, end, mc, cov]
            ll[1] = str(int(ll[1]) - 1)  # because bed is 0 based
            bed_line = '\t'.join([chrom, ll[1], ll[1], ll[4], ll[5]]) + '\n'
            # assign each line to its patten content,
            # will write multiple times if patten overlap
            for c, p in pattern_dict.items():
                if ll[3] in p:
                    handle_dict[c].write(bed_line)
    # close handle
    for handle in handle_dict.values():
        handle.close()
    return path_dict


def map_to_region(allc_path, out_path_prefix,
                  region_bed_path, region_name, genome_size_path,
                  context_pattern, max_cov_cutoff, remove_tmp):
    """
    Map one allc file into many region set bed file using bedtools map.
    Count mC and coverage in each region for each context pattern.

    Parameters
    ----------
    allc_path
    out_path_prefix
    region_bed_path
    region_name
    genome_size_path
        UCSC chrom.sizes file, will determine which chrom to keep in the output.
        Use main chrom if want to remove those random contigs
    context_pattern
    max_cov_cutoff
    remove_tmp

    Returns
    -------

    """

    # parse ref chrom with ordered chromosome
    ref_chrom_dict = parse_chrom_size(genome_size_path)

    # prepare ALLC bed dict, split ALLC into different contexts
    # bed format [chrom, start, end, mc, cov]
    print('Splitting ALLC')
    # split chromosome and avoid sorting
    allc_bed_path_dict = _split_to_chrom_bed(allc_path=allc_path,
                                             context_pattern=context_pattern,
                                             out_path_prefix=out_path_prefix + '.tmp',
                                             genome_size_path=genome_size_path,
                                             max_cov_cutoff=max_cov_cutoff)
    # concat bed with ordered chromosome
    tmp_dict = {}
    for c in context_pattern:
        c_path_list = [allc_bed_path_dict[(c, _chrom)]
                       for _chrom in ref_chrom_dict.keys()
                       if (c, _chrom) in allc_bed_path_dict]
        cmd = ['cat'] + c_path_list
        concat_bed_path = out_path_prefix + f'.{c}.tmp.total.bed'
        with open(concat_bed_path, 'w') as fh:
            run(cmd, stdout=fh)
            for p in c_path_list:
                run(['rm', '-f', p])
        tmp_dict[c] = concat_bed_path
    allc_bed_path_dict = tmp_dict

    print('Reading ALLC Bed')
    allc_bed_dict = {k: BedTool(path) for k, path in allc_bed_path_dict.items()}
    # k is (context_pattern)

    # prepare all region bed files
    print('Reading Region Bed')
    if len(region_bed_path) != len(region_name):
        raise ValueError('Number of region BED path != Number of region names')
    # input region bed, sort across allc chrom order
    # chrom_order is in UCSC genome size format,
    # make a bed format from it and then intersect with bed_p to filter out chromosomes not appear in ALLC

    region_bed_dict = {region_n: BedTool(bed_p).sort(g=genome_size_path)
                       for bed_p, region_n in zip(region_bed_path, region_name)}

    # bedtools map
    for context_name, allc_bed in allc_bed_dict.items():
        for region_name, region_bed in region_bed_dict.items():
            print(f'Map {context_name} ALLC Bed to {region_name} Region Bed')
            region_bed.map(b=allc_bed, c='4,5', o='sum,sum', g=genome_size_path) \
                .saveas(out_path_prefix + f'.{region_name}_{context_name}.count_table.bed.gz',
                        compressed=True)

    # cleanup the tmp bed files.
    if remove_tmp:
        print('Clean tmp Bed file')
        for path in allc_bed_path_dict.values():
            run(['rm', '-f', path])
    cleanup()  # pybedtools tmp files
    print('Finish')
    return


def merge_allc(allc_paths, out_path, cpu=1, index=False,
               get_mcg=True, cg_pattern='CGN'):
    """
    Just a wrapper of methylpy merge allc
    :param allc_paths:
    :param out_path:
    :param cpu:
    :param index:
    :param get_mcg:
    :param cg_pattern:
    :return:
    """
    if len(allc_paths) == 1:
        with open(allc_paths[0]) as f:
            allc_paths = [i.strip('\n') for i in f.readlines()]
    if not os.path.exists(out_path + '.gz'):
        merge_allc_files(allc_paths,
                         out_path,
                         num_procs=cpu,
                         mini_batch=150,
                         compress_output=False,
                         skip_snp_info=True,
                         buffer_line_number=100000,
                         index=False)
        # use bgzip and tabix
        run(['bgzip', out_path])

        if index:
            run(['tabix', '-b', '2', '-e', '2', '-s', '1', out_path + '.gz'])
    if get_mcg:
        extract_mcg(allc_path=out_path + '.gz', out_path=out_path[:-6] + f'{cg_pattern}.tsv.gz', cg_pattern=cg_pattern)
    return


def allc_to_bigwig(allc_path, out_path, chrom_size, mc_type='CGN'):
    from .methylpy_utilities import convert_allc_to_bigwig
    convert_allc_to_bigwig(allc_path,
                           out_path,
                           chrom_size,
                           mc_type=mc_type,
                           bin_size=100,
                           path_to_wigToBigWig="",
                           path_to_samtools="",
                           min_bin_sites=0,
                           min_bin_cov=0,
                           max_site_cov=None,
                           min_site_cov=0,
                           add_chr_prefix=True
                           )
    return


def extract_mcg(allc_path, out_path, merge_strand=True, header=False, cg_pattern='CGN'):
    if '.gz' in allc_path[-3:]:
        opener = partial(gzip.open, mode='rt')
    else:
        opener = partial(open, mode='r')
    writer = partial(gzip.open, mode='wt')

    context_set = parse_mc_pattern(cg_pattern)
    with opener(allc_path) as allc, \
            writer(out_path) as out_allc:
        if header:
            allc.readline()
        if merge_strand:
            prev_line = None
            cur_chrom = None
            for line in allc:
                cur_line = line.strip('\n').split('\t')
                if cur_line[3] not in context_set:
                    continue
                if cur_line[0] != cur_chrom:
                    if prev_line is not None:
                        out_allc.write('\t'.join(prev_line) + '\n')
                    prev_line = cur_line
                    cur_chrom = cur_line[0]
                    continue
                if prev_line is None:
                    prev_line = cur_line
                    continue
                else:
                    # pos should be continuous, strand should be reverse
                    if int(prev_line[1]) + 1 == int(cur_line[1]) and prev_line[2] != cur_line[2]:
                        new_line = prev_line[:4] + [str(int(prev_line[4]) + int(cur_line[4])),
                                                    str(int(prev_line[5]) + int(cur_line[5])), '1']
                        out_allc.write('\t'.join(new_line) + '\n')
                        prev_line = None
                    # otherwise, only write and update prev_line
                    else:
                        out_allc.write('\t'.join(prev_line) + '\n')
                        prev_line = cur_line
        else:
            for line in allc:
                cur_line = line.strip('\n').split('\t')
                if cur_line[3] not in context_set:
                    continue
                out_allc.write('\t'.join(cur_line) + '\n')
    print('Extract CG finished:', out_path)
    return


def get_allc_profile(allc_path, drop_n=True, n_rows=100000000, out_path=None):
    """
    Generate approximate profile for allc file. 1e8 rows finish in about 5 min.

    Parameters
    ----------
    allc_path
        path of the allc file
    drop_n
        whether drop context contain N
    n_rows
        number of rows to use, 1e8 is sufficient to get an approximate profile
    out_path
        if not None, save profile to out_path
    Returns
    -------

    """
    if 'gz' in allc_path:
        opener = partial(gzip.open, mode='rt')
    else:
        opener = partial(open, mode='r')

    # initialize count dict
    mc_sum_dict = defaultdict(int)
    cov_sum_dict = defaultdict(int)
    cov_sum2_dict = defaultdict(int)  # sum of square, for calculating variance
    rate_sum_dict = defaultdict(float)
    rate_sum2_dict = defaultdict(float)  # sum of square, for calculating variance
    context_count_dict = defaultdict(int)
    with opener(allc_path) as f:
        n = 0
        for line in f:
            chrom, pos, strand, context, mc, cov, p = line.split('\t')
            if drop_n and 'N' in context:
                continue
            # mc and cov
            mc_sum_dict[context] += int(mc)
            cov_sum_dict[context] += int(cov)
            cov_sum2_dict[context] += int(cov) ** 2
            # raw base rate
            rate = int(mc) / int(cov)
            rate_sum_dict[context] += rate
            rate_sum2_dict[context] += rate ** 2
            # count context finally
            context_count_dict[context] += 1
            n += 1
            if (n_rows is not None) and (n >= n_rows):
                break
    # overall count
    profile_df = pd.DataFrame({'partial_mc': mc_sum_dict,
                               'partial_cov': cov_sum_dict})
    profile_df['base_count'] = pd.Series(context_count_dict)
    profile_df['overall_mc_rate'] = profile_df['partial_mc'] / profile_df['partial_cov']

    # cov base mean and base std.
    # assume that base cov follows normal distribution
    cov_sum_series = pd.Series(cov_sum_dict)
    cov_sum2_series = pd.Series(cov_sum2_dict)
    profile_df['base_cov_mean'] = cov_sum_series / profile_df['base_count']
    profile_df['base_cov_std'] = np.sqrt(
        (cov_sum2_series / profile_df['base_count']) - profile_df['base_cov_mean'] ** 2)

    # assume that base rate follow beta distribution
    # so that observed rate actually follow joint distribution of beta (rate) and normal (cov) distribution
    # here we use the observed base_rate_mean and base_rate_var to calculate
    # approximate alpha and beta value for the base rate beta distribution
    rate_sum_series = pd.Series(rate_sum_dict)
    rate_sum2_series = pd.Series(rate_sum2_dict)
    profile_df['base_rate_mean'] = rate_sum_series / profile_df['base_count']
    profile_df['base_rate_var'] = (rate_sum2_series / profile_df['base_count']) - profile_df['base_rate_mean'] ** 2

    # based on beta distribution mean, var
    # a / (a + b) = base_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = base_rate_var
    # we have:
    a = (1 - profile_df['base_rate_mean']) * (profile_df['base_rate_mean'] ** 2) / profile_df['base_rate_var'] - \
        profile_df['base_rate_mean']
    b = a * (1 / profile_df['base_rate_mean'] - 1)
    profile_df['base_beta_a'] = a
    profile_df['base_beta_b'] = b

    if out_path is not None:
        profile_df.to_csv(out_path, sep='\t')
        return None
    else:
        return profile_df
