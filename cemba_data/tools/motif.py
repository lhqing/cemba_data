import subprocess
import collections
import shlex
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from .utilities import parse_chrom_size, parse_mc_pattern
from .allc._open import open_allc, open_gz
import logging

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def get_motif_bed(motif_file_path, genome, out_prefix, non_overlap=True, non_overlap_kws=None,
                  score_cutoff=8, cpu=1, keep_all=True, homer_path=None):
    """
    Get a motif bed file [chrom, start, end] from homer motif matrix file, using homer scanMotifGenomeWide.pl
    The output is bgzip & tabixed, and usually non-overlapped,
    because get_motif_rate must have non-overlap input.
    See homer doc here: http://homer.ucsd.edu/homer/motif/genomeWideMotifScan.html

    Parameters
    ----------
    motif_file_path
        Path to motif matrix file, homer format.
        See doc here: http://homer.ucsd.edu/homer/motif/creatingCustomMotifs.html
    genome
        this directly goes to homer, accept genome name like 'mm10' or a fasta file path
    out_prefix
        prefix of the output, suffix will be added automatically
    non_overlap
        If true, will use bedtools merge to merge all overlap regions
    non_overlap_kws
        dict key word parameters pass to bedtools merge
    score_cutoff
        homer motif match score used to filter motif regions
    cpu
        Number of CPU to use in homer
    keep_all
        Homer parameter, keep ALL sites, even ones that overlap
    homer_path
        Path to homer executable files
    Returns
    -------

    """
    if homer_path is not None:
        os.environ['PATH'] += ':' + homer_path
    # check if homer function can be found
    try:
        subprocess.run(['scanMotifGenomeWide.pl'], check=True)
    except FileNotFoundError:
        raise FileNotFoundError('Homer file scanMotifGenomeWide.pl can not found, '
                                'add homer to PYTHONPATH or provide '
                                '/path/to/homer/bin/ in homer_path parameter')

    _keep_all = ' '
    if keep_all:
        _keep_all = ' -keepAll'
    cmd = f'scanMotifGenomeWide.pl {motif_file_path} {genome}{_keep_all} -bed -p {cpu} -int'

    score_count = collections.defaultdict(int)
    total_line_count = 0
    pass_line_count = 0
    out_path = out_prefix + '.motif.bed'
    with open(out_path, 'w') as f:
        process = subprocess.Popen(shlex.split(cmd),
                                   encoding='utf8',
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        for line in process.stdout:
            # the homer original output example:
            # chr10   3101568 3101577 Lhx1(Homeobox)  7       +
            # chr10   3101570 3101579 Lhx1(Homeobox)  7       -
            # chr10   3103413 3103422 Lhx1(Homeobox)  6       +
            # chr10   3103415 3103424 Lhx1(Homeobox)  6       -
            # This function only keep chrom, start, end.
            ll = line.strip().split('\t')
            score = int(ll[4])
            score_count[score] += 1
            total_line_count += 1
            if score >= score_cutoff:
                new_line = '\t'.join(ll[:3]) + '\n'
                f.write(new_line)
                pass_line_count += 1

    print(f'Found a total of {total_line_count} motifs, {pass_line_count} have score >= {score_cutoff}.')
    score_count = collections.OrderedDict(sorted(score_count.items(), key=lambda t: t[0]))
    for k, v in score_count.items():
        print(f'Score: {k}; Count: {v}')

    if non_overlap:
        if non_overlap_kws is None:
            non_overlap_kws = {}
        import pybedtools
        non_overlap_bed = pybedtools.BedTool(out_path).merge(**non_overlap_kws).saveas(out_path)

        region_count = non_overlap_bed.count()
        print(f'Make bed file intervals Non-overlap, after bedtools merge, '
              f'{region_count} regions remained.')

    p = subprocess.run(['bgzip', out_path], stderr=subprocess.PIPE, encoding='utf8')
    if p.returncode != 0:
        raise OSError(p.stderr)
    p = subprocess.run(['tabix', out_path + '.gz'], stderr=subprocess.PIPE, encoding='utf8')
    if p.returncode != 0:
        raise OSError(p.stderr)


def _check_bgz_bed_non_overlap(bed_path):
    if not bed_path.endswith('.gz'):
        raise ValueError('Bed file should be bgzip/tabixed to use this function.')

    cur_chrom = 'NOTACHROM'
    last_end = -1
    line_number = 0
    with open_gz(bed_path) as f:
        for line in f:
            line_number += 1
            chrom, start, end, *_ = line.strip('\n').split('\t')
            if chrom == cur_chrom:
                start = int(start)
                end = int(end)
                if start >= last_end:
                    last_end = end
                else:
                    return False
            else:
                cur_chrom = chrom
                last_end = -1
    return True


def _get_tabix_chroms(file_path):
    process = subprocess.run(['tabix', file_path, '-l'], check=True,
                             stdout=subprocess.PIPE, encoding='utf8')

    chroms = process.stdout.strip().split('\n')
    return chroms


def _get_non_overlap_motif_rate(allc_path, bed_path, chrom_size_path, mc_types):
    mc_sets = {mc_type: parse_mc_pattern(mc_type) for mc_type in mc_types}

    if not os.path.exists(allc_path + '.tbi'):
        raise FileNotFoundError(f'ALLC file {allc_path} do not have .tbi index')
    if not os.path.exists(bed_path + '.tbi'):
        raise FileNotFoundError(f'Bed file {bed_path} do not have .tbi index')

    ref_chrom_dict = parse_chrom_size(chrom_size_path)
    allc_chroms = _get_tabix_chroms(allc_path)
    bed_chroms = _get_tabix_chroms(bed_path)
    use_chroms = (set(ref_chrom_dict.keys()) & set(allc_chroms) & set(bed_chroms))

    mc_count_dict = {mc_type: 0 for mc_type in mc_types}
    cov_count_dict = {mc_type: 0 for mc_type in mc_types}
    for chrom in sorted(list(use_chroms)):
        # loop on each chrom
        with open_allc(allc_path, region=chrom) as allc, \
                open_gz(bed_path, region=chrom) as bed:
            start = -1
            end = -1
            for line in allc:
                _, pos, strand, context, mc, cov, _ = line.strip().split('\t')
                pos = int(pos)
                if pos < start:
                    # not reach the region
                    continue
                elif pos <= end:  # note: bed format is right close
                    # in the region
                    for mc_type, mc_set in mc_sets.items():
                        if context in mc_set:
                            mc_count_dict[mc_type] += int(mc)
                            cov_count_dict[mc_type] += int(cov)
                else:
                    # out of the region
                    while pos > end:
                        # read new regions until the end >= current pos
                        try:
                            _, start, end = bed.readline().strip().split('\t')
                        except ValueError:
                            # read to the end of the file,
                            # no need to check remaining lines
                            break
                        start = int(start)
                        end = int(end)
                    if pos > start:
                        # check this line and new region
                        for mc_type, mc_set in mc_sets.items():
                            if context in mc_set:
                                mc_count_dict[mc_type] += int(mc)
                                cov_count_dict[mc_type] += int(cov)
    records = []
    for mc_type in mc_types:
        mc_count = mc_count_dict[mc_type]
        cov_count = cov_count_dict[mc_type]
        records.append({'mc_type': mc_type,
                        'mc': mc_count, 'cov': cov_count})
    return pd.DataFrame(records)


def batch_get_motif_rate(allc_path_table, bed_path, chrom_size_path,
                         out_path, mc_types, cpu=10):
    """
    Map a list of ALLC files into a non-overlap bed file.
    But only count the overall mC and cov across the whole bed file.

    Parameters
    ----------
    allc_path_table
        Two column table without header. The 1st is file id, the 2nd is the allc path
    bed_path
        Path of the non-overlap bed file
    chrom_size_path
        Path of the UCSC chrom size file
    out_path
        Path of the output table
    mc_types
        mC context pattern(s)
    cpu
        Number of CPU to use

    Returns
    -------

    """
    if not _check_bgz_bed_non_overlap(bed_path):
        raise ValueError('Bed file is not non-overlap, this function only take non-overlap intervals. '
                         'Use bedtools merge to merge that. Or set non_overlap=True in get_motif_bed function.')
    cell_allc_paths = pd.read_table(allc_path_table, index_col=0, squeeze=True, header=None)
    if not isinstance(cell_allc_paths, pd.Series):
        raise ValueError('allc_path_table is malformed. It should only contain 2 columns without header line, '
                         'the 1st is file id, the 2nd is the allc path.')

    if isinstance(mc_types, str):
        mc_types = [mc_types]

    with ProcessPoolExecutor(max_workers=cpu) as executor:
        future_result = {executor.submit(_get_non_overlap_motif_rate,
                                         allc_path=allc_path,
                                         bed_path=bed_path,
                                         chrom_size_path=chrom_size_path,
                                         mc_types=mc_types): file_id
                         for file_id, allc_path in cell_allc_paths.iteritems()}
        records = []
        for future in as_completed(future_result):
            file_id = future_result[future]
            try:
                result = future.result()
                result['file_id'] = file_id
                records.append(result)
            except Exception as exc:
                log.info('%r generated an exception: %s' % (file_id, exc))
    total_result = pd.concat(records, sort=True, ignore_index=True)
    total_result.to_csv(out_path, sep='\t', index=None)
    return
