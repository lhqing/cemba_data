from .utilities import *
import gzip
from functools import partial
from pybedtools import BedTool, cleanup
from pybedtools.helpers import BEDToolsError
from subprocess import run
import argparse
import os
from pandas import Series


def split_to_bed(allc_path, context_pattern, genome_size_path,
                 out_path_prefix, max_cov_cutoff=None,
                 split_chromosome=False,
                 compression=True, gzip_level=3):
    """
    Split ALLC into bed format, chrom column contain "chr".
    :param allc_path: Single ALLC file path
    :param context_pattern: comma separate context patterns or list
    :param out_path_prefix: Single output prefix
    :param max_cov_cutoff: 2 for single cell, None for bulk or merged allc
    :param split_chromosome: 
    :param genome_size_path: UCSC chrom size file path
    :param compression: gzip compression or not
    :param gzip_level: compression level, default 3
    :return: Path dict for out put files
    """
    # check whether chr is in chrom:
    with gzip.open(allc_path, 'rt') as allc:
        test_line = allc.readline()
    add_chr = 'chr' != test_line[:3]
    cur_chrom = test_line.split('\t')[0]
    if add_chr:
        cur_chrom = 'chr' + cur_chrom
    chrom_order_list = [cur_chrom]

    # prepare context
    if isinstance(context_pattern, str):
        context_pattern = context_pattern.split(',')
    pattern_dict = {c: parse_mc_pattern(c) for c in context_pattern}

    # prepare out path
    if split_chromosome:
        chrom_list = list(parse_chrom_szie(genome_size_path).keys())
        path_dict = {(c, chrom): out_path_prefix + f'.{c}.{chrom}.bed'
                     for c in context_pattern
                     for chrom in chrom_list}
    else:
        path_dict = {c: out_path_prefix + f'.{c}.bed' for c in context_pattern}

    # open all paths
    if compression:
        open_func = partial(gzip.open, mode='wt', compresslevel=gzip_level)
        path_dict = {k: v + '.gz' for k, v in path_dict.items()}  # add gz to path
    else:
        open_func = partial(open, mode='w')

    handle_dict = {k: open_func(path) for k, path in path_dict.items()}

    # split ALLC
    with gzip.open(allc_path, 'rt') as allc:
        for line in allc:
            ll = line.split('\t')
            # filter max cov (for single cell data)
            if max_cov_cutoff is not None:
                if int(ll[5]) > max_cov_cutoff:
                    continue
            # bed format [chrom, start, end, mc, cov]
            ll[1] = str(int(ll[1]) - 1)  # because bed is 0 based
            bed_line = '\t'.join([ll[0], ll[1], ll[1], ll[4], ll[5]]) + '\n'
            # add "chr" to chrom
            if add_chr:
                bed_line = 'chr' + bed_line
                chrom = 'chr' + ll[0]
            else:
                chrom = ll[0]
            if chrom == 'chrL':
                continue
            # record ALLC chrom order
            if chrom != chrom_order_list[-1]:
                chrom_order_list.append(chrom)
            # assign each line to its patten content
            for c, p in pattern_dict.items():
                if ll[3] in p:
                    if split_chromosome:
                        handle_dict[(c, chrom)].write(bed_line)
                    else:
                        handle_dict[c].write(bed_line)

    # close handle
    for handle in handle_dict.values():
        handle.close()
    with open(out_path_prefix + '.chrom_order', 'w') as f:
        chrom_series = Series(parse_chrom_szie(genome_size_path)).reindex(chrom_order_list).dropna()
        for chrom, length in chrom_series.iteritems():
            f.write(f'{chrom}\t{length}\n')
    return path_dict


def map_to_region(allc_path, out_path_prefix,
                  region_bed_path, region_name, genome_size_path,
                  context_pattern, max_cov_cutoff,
                  remove_tmp=True, tmp_compression=False):
    """
    Map one allc file into many region set bed file using bedtools map.
    Count mC and coverage in each region for each context pattern.
    :param allc_path:
    :param out_path_prefix:
    :param region_bed_path:
    :param region_name:
    :param genome_size_path:
    :param context_pattern:
    :param max_cov_cutoff:
    :param remove_tmp:
    :param tmp_compression:
    :return:
    """
    # prepare ALLC bed dict, split ALLC into different contexts
    # bed format [chrom, start, end, mc, cov]
    print('Splitting ALLC')
    # split chromosome and avoid sorting
    allc_bed_path_dict = split_to_bed(allc_path=allc_path,
                                      context_pattern=context_pattern,
                                      out_path_prefix=out_path_prefix + '.tmp',
                                      split_chromosome=False,
                                      genome_size_path=genome_size_path,
                                      max_cov_cutoff=max_cov_cutoff,
                                      compression=tmp_compression)
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
    with open(out_path_prefix + '.tmp.chrom_order', 'r') as chrom_size, \
            open(out_path_prefix + '.tmp.chrom_order.bed', 'w') as chrom_bed:
        for line in chrom_size:
            ll = line.strip('\n').split('\t')
            chrom_bed.write('\t'.join([ll[0], '0', ll[1]]) + '\n')
    allc_chrom_bed = BedTool(out_path_prefix + '.tmp.chrom_order.bed')  # bed file for allc chrom
    region_bed_dict = {region_n: BedTool(bed_p).intersect(allc_chrom_bed).sort(g=out_path_prefix + '.tmp.chrom_order')
                       for bed_p, region_n in zip(region_bed_path, region_name)}

    # bedtools map
    for context_name, allc_bed in allc_bed_dict.items():
        for region_name, region_bed in region_bed_dict.items():
            print(f'Map {context_name} ALLC Bed to {region_name} Region Bed')
            try:
                region_bed.map(b=allc_bed, c='4,5', o='sum,sum', g=out_path_prefix + '.tmp.chrom_order') \
                    .sort(g=genome_size_path) \
                    .saveas(out_path_prefix + f'.{region_name}_{context_name}.count_table.bed.gz',
                            compressed=True)
            except BEDToolsError:
                # for some rare cases (10 in 2700),
                # the gzip bed raise error in bedtools map, but unzip bed works fine.
                # https://github.com/arq5x/bedtools2/issues/363
                # not sure why, but try unzip once
                if tmp_compression:
                    print('Encounter BEDToolsError')
                    if os.path.exists(allc_bed_path_dict[context_name]):
                        # unzip the fine if it haven't
                        run(['gunzip', allc_bed_path_dict[context_name]])
                    allc_bed = BedTool(allc_bed_path_dict[context_name][:-3])  # open unzip bed
                    region_bed.map(b=allc_bed, c='4,5', o='sum,sum', g=out_path_prefix + '.tmp.chrom_order') \
                        .sort(g=genome_size_path) \
                        .saveas(out_path_prefix + f'.{region_name}_{context_name}.count_table.bed.gz',
                                compressed=True)
                    print('Solved BEDToolsError when gunzip the file...')
                else:
                    # file already not zipped, some novel error...
                    raise BEDToolsError

    # cleanup the tmp bed files.
    if remove_tmp:
        print('Clean tmp Bed file')
        for path in allc_bed_path_dict.values():
            if os.path.exists(path):
                run(['rm', '-f', path])
            else:  # bed have been unzipped
                run(['rm', '-f', path[:-3]])
        run(['rm', '-f', out_path_prefix + '.tmp.chrom_order'])
        run(['rm', '-f', out_path_prefix + '.tmp.chrom_order.bed'])

    cleanup()  # pybedtools tmp files
    print('Finish')
    return


def map_to_region_register_subparser(subparser):
    parser = subparser.add_parser('map-to-region',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Map ALLC file to region BED, "
                                       "get base mC and coverage count for each region.")
    parser.set_defaults(func=map_to_region)

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="One ALLC file path"
    )

    parser_req.add_argument(
        "--out_path_prefix",
        type=str,
        required=True,
        help="Output file prefix"
    )

    parser_req.add_argument(
        "--region_bed_path",
        type=str,
        required=True,
        nargs='+',
        help="Space separated region BED file paths"
    )

    parser_req.add_argument(
        "--region_name",
        type=str,
        required=True,
        nargs='+',
        help="Space separated region set names corresponding to --region_bed_path"
    )

    parser_req.add_argument(
        "--genome_size_path",
        type=str,
        required=True,
        help="UCSC genome size file"
    )

    parser_req.add_argument(
        "--context_pattern",
        type=str,
        required=True,
        nargs='+',
        help="Space separated methylation context pattern, N for ATCG, H for ATC"
    )

    parser_opt.add_argument(
        "--max_cov_cutoff",
        type=int,
        required=False,
        default=None,
        help="Maximum cutoff for coverage in each base, "
             "e.g. 2 for single cell data, None for bulk seq."
    )

    parser_opt.add_argument(
        "--remove_tmp",
        type=bool,
        required=False,
        default=True,
        help="Remove tmp file or not"
    )

    parser_opt.add_argument(
        "--tmp_compression",
        type=bool,
        required=False,
        default=False,
        help="Compress tmp file (slower but space efficient) or not"
    )
    return
