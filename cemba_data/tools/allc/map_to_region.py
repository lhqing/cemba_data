from typing import Union, Tuple
from collections import defaultdict
import subprocess
import shlex
from functools import partial
from .utilities import check_tbi_chroms
from ._open import open_allc, open_gz
from ..utilities import parse_chrom_size, parse_mc_pattern


def _allc_to_site_bed(allc_path: str,
                      out_prefix: str,
                      chrom_size_path: str,
                      mc_contexts: Union[str, list],
                      max_cov_cutoff: int = 9999) -> Tuple[list, list]:
    """
    1. Split allc by context into several site-bed files
    2. Keep chromosome sorted as the chrom_size_path
    3. BED file 4th column is mc, 5th column is cov, and each row is a single site
    """
    # TODO add parallel to this
    out_prefix = out_prefix.rstrip('.')
    if isinstance(mc_contexts, str):
        mc_contexts = mc_contexts.split(' ')
    mc_contexts = list(set(mc_contexts))

    # because mc_contexts can overlap (e.g. CHN, CAN)
    # each context may associate to multiple handle
    context_handle = defaultdict(list)
    handle_collect = []
    out_paths = []
    for mc_context in mc_contexts:
        out_path = out_prefix + f'.extract_{mc_context}.bed.gz'
        out_paths.append(out_path)
        _handle = open_allc(out_path, 'w')
        handle_collect.append(_handle)
        parsed_context_set = parse_mc_pattern(mc_context)
        for pattern in parsed_context_set:
            context_handle[pattern].append(_handle)

    # split file first
    chrom_size_dict = parse_chrom_size(chrom_size_path)
    with open_allc(allc_path, region=' '.join(chrom_size_dict.keys())) as allc:
        for line in allc:
            chrom, pos, _, context, mc, cov, *_ = line.strip('\n').split('\t')
            if int(cov) > max_cov_cutoff:
                continue
            bed_line = '\t'.join([chrom, pos, pos, mc, cov]) + '\n'
            try:
                [h.write(bed_line) for h in context_handle[context]]
            except KeyError:
                continue
    for handle in handle_collect:
        handle.close()
    return mc_contexts, out_paths


def _bedtools_map(region_bed, site_bed, out_bed, save_zero_cov=True):
    cmd = f'bedtools map -a {region_bed} -b {site_bed} -c 4,5 -o sum,sum'
    bed_out = subprocess.run(shlex.split(cmd),
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             encoding='utf8', check=True)
    if out_bed.endswith('gz'):
        opener = partial(open_gz, mode='wt')
    else:
        opener = partial(open, mode='w')

    with opener(out_bed) as out_handle:
        for line in bed_out.stdout:
            if save_zero_cov:
                out_handle.write(line)
            else:
                # the last item is cov
                if not line.endswith('\t0\n'):
                    out_handle.write(line)
    return


def allc_to_region_count(allc_path,
                         out_prefix,
                         region_bed_paths,
                         region_bed_names,
                         mc_contexts,
                         chrom_size_path,
                         max_cov_cutoff,
                         save_zero_cov=True,
                         remove_tmp=True):
    # 1. bgzip
    # 2. order of chrom should be the same as order of chrom_size_path
    genome_dict = parse_chrom_size(chrom_size_path)
    for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
        if not check_tbi_chroms(region_bed_path, genome_dict):
            raise ValueError(f'The bed file {region_bed_path} chromosome order is different '
                             f'from the {chrom_size_path}')

    print('Extract ALLC context')
    out_prefix = out_prefix.rstrip('.')
    mc_contexts, site_bed_paths = _allc_to_site_bed(allc_path=allc_path,
                                                    out_prefix=out_prefix,
                                                    chrom_size_path=chrom_size_path,
                                                    mc_contexts=mc_contexts,
                                                    max_cov_cutoff=max_cov_cutoff)
    print('Map to regions')
    save_flag = 'full' if save_zero_cov else 'sparse'
    for region_name, region_bed_path in zip(region_bed_names, region_bed_paths):
        for mc_context, site_bed_path in zip(mc_contexts, site_bed_paths):
            try:
                _bedtools_map(region_bed=region_bed_path,
                              site_bed=site_bed_path,
                              out_bed=out_prefix+f'.{region_name}_{mc_context}.{save_flag}.bed.gz',
                              save_zero_cov=save_zero_cov)
            except subprocess.CalledProcessError as e:
                print(e.stderr)
                raise e

    if remove_tmp:
        for site_bed_path in site_bed_paths:
            subprocess.run(['rm', '-f', site_bed_path])
    return
