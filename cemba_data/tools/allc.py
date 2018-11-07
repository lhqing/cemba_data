from .utilities import *
import gzip
from functools import partial
from pybedtools import BedTool, cleanup
from pybedtools.helpers import BEDToolsError
from subprocess import run
import argparse
import os
from .methylpy_utilities import merge_allc_files


def _split_to_chrom_bed(allc_path, context_pattern, genome_size_path,
                        out_path_prefix, max_cov_cutoff=None,
                        compression=False, gzip_level=2, remove_chrm=True):
    """
    Split ALLC into bed format, chrom column contain "chr".
    :param allc_path: Single ALLC file path
    :param context_pattern: comma separate context patterns or list
    :param out_path_prefix: Single output prefix
    :param max_cov_cutoff: 2 for single cell, None for bulk or merged allc
    :param genome_size_path: UCSC chrom size file path
    :param compression: gzip compression or not
    :param gzip_level: compression level, default 3
    :return: Path dict for out put files
    """
    if remove_chrm:
        remove_chrm = ['chrM']
    else:
        remove_chrm = None
    chrom_set = set(parse_chrom_size(genome_size_path,
                                     remove_chr_list=remove_chrm,
                                     add_chr=add_chr).keys())
    # prepare context
    if isinstance(context_pattern, str):
        context_pattern = context_pattern.split(',')
    pattern_dict = {c: parse_mc_pattern(c) for c in context_pattern}

    # prepare out path
    path_dict = {(c, chrom): out_path_prefix + f'.{c}.{chrom}.bed'
                 for c in context_pattern
                 for chrom in chrom_set}
    # open all paths
    if compression:
        open_func = partial(gzip.open, mode='at', compresslevel=gzip_level)
        path_dict = {k: v + '.gz' for k, v in path_dict.items()}  # add gz to path
    else:
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
    _add_chr = False
    with allc_open_func(allc_path) as allc:
        for line in allc:
            if first:
                chrom = line.split('\t')[0]
                _add_chr = add_chr & ('chr' != line[:3])
                if _add_chr:
                    chrom = 'chr' + chrom
                if chrom not in chrom_set:
                    continue
                first = False
                cur_chrom = chrom
                handle_dict = {c: open_func(path_dict[(c, cur_chrom)]) for c in context_pattern}
            ll = line.split('\t')
            # filter max cov (for single cell data)
            if max_cov_cutoff is not None:
                if int(ll[5]) > max_cov_cutoff:
                    continue
            # add "chr" to chrom
            if _add_chr:
                chrom = 'chr' + ll[0]
            else:
                chrom = ll[0]
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
                  context_pattern, max_cov_cutoff,
                  remove_tmp, tmp_compression, add_chr):
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
    :param add_chr:
    :return:
    """
    # parse ref chrom with ordered chromosome
    ref_chrom_dict = parse_chrom_size(genome_size_path, add_chr=add_chr)

    # prepare ALLC bed dict, split ALLC into different contexts
    # bed format [chrom, start, end, mc, cov]
    print('Splitting ALLC')
    # split chromosome and avoid sorting

    allc_bed_path_dict = _split_to_chrom_bed(allc_path=allc_path,
                                             context_pattern=context_pattern,
                                             out_path_prefix=out_path_prefix + '.tmp',
                                             genome_size_path=genome_size_path,
                                             max_cov_cutoff=max_cov_cutoff,
                                             compression=tmp_compression,
                                             add_chr=add_chr)
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
            try:
                region_bed.map(b=allc_bed, c='4,5', o='sum,sum', g=genome_size_path) \
                    .saveas(out_path_prefix + f'.{region_name}_{context_name}.count_table.bed.gz',
                            compressed=True)
            except BEDToolsError as err:
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
                    region_bed.map(b=allc_bed, c='4,5', o='sum,sum', g=genome_size_path) \
                        .saveas(out_path_prefix + f'.{region_name}_{context_name}.count_table.bed.gz',
                                compressed=True)
                    print('Solved BEDToolsError when gunzip the file...')
                else:
                    # file already unzipped, some novel error...
                    raise err

    # cleanup the tmp bed files.
    if remove_tmp:
        print('Clean tmp Bed file')
        for path in allc_bed_path_dict.values():
            if os.path.exists(path):
                run(['rm', '-f', path])
            else:  # bed have been unzipped
                run(['rm', '-f', path[:-3]])

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

    parser_opt.add_argument(
        "--add_chr",
        type=bool,
        required=False,
        default=False,
        help="add 'chr' before chromosome or not"
    )
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


def allc_to_bigwig_register_subparser(subparser):
    parser = subparser.add_parser('allc-to-bigwig',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Just a wrapper of methylpy allc-to-bigwig")
    parser.set_defaults(func=allc_to_bigwig)

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_path",
        type=str,
        required=True,
        help="ALLC path"
    )
    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Out path"
    )
    parser_req.add_argument(
        "--chrom_size",
        type=str,
        required=True,
        help="UCSC chrom.sizes format indicating genome size. ALLC Chr not in this file will be removed."
    )

    parser_opt.add_argument(
        "--mc_type",
        type=str,
        required=False,
        default='CGN',
        help="mC context pattern to use"
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


def merge_allc_register_subparser(subparser):
    parser = subparser.add_parser('merge-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Just a wrapper of methylpy merge_allc_files, "
                                       "without doing methylpy's index. "
                                       "But use bgzip and tabix instead")
    parser.set_defaults(func=merge_allc)

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--allc_paths",
        type=str,
        required=True,
        nargs='+',
        help="Space separated ALLC paths OR a file contains all ALLC paths rows. "
             "If provide only one path, each row in the file should be a path"
    )

    parser_req.add_argument(
        "--out_path",
        type=str,
        required=True,
        help="Output path for the merged ALLC file"
    )

    parser_opt.add_argument(
        "--cpu",
        type=int,
        required=False,
        default=1,
        help="Number of CPUs for merge ALLC, parallel on chrom bins level."
    )

    parser_opt.add_argument(
        "--index",
        type=bool,
        required=False,
        default=False,
        help="methylpy default index, not doing it by default."
    )

    parser_opt.add_argument(
        "--get_mcg",
        type=bool,
        required=False,
        default=True,
        help="Generate a CG only, strand merged ALLC file for CG DMR calling."
    )

    parser_opt.add_argument(
        "--cg_pattern",
        type=str,
        required=False,
        default='CGN',
        help="mCG context pattern for --get_mcg."
    )
