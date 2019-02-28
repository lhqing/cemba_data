import subprocess
import collections
import shlex
import os
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from .utilities import parse_chrom_size, parse_mc_pattern
from .open import open_allc, open_gz
import logging

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def get_motif_bed(motif_file_path, genome, out_prefix, non_overlap=True, non_overlap_kws=None,
                  score_cutoff=8, bed_out=True, cpu=1, keep_all=True, homer_path=None):
    if homer_path is not None:
        os.environ['PATH'] += ':' + homer_path
    # check if homer function can be found
    try:
        subprocess.run(['scanMotifGenomeWide.pl'], check=True)
    except FileNotFoundError:
        raise FileNotFoundError('Homer file scanMotifGenomeWide.pl can not found, '
                                'add homer to PYTHONPATH or provide '
                                '/path/to/homer/bin/ in homer_path parameter')

    _bed_out = ' '
    if bed_out:
        _bed_out = ' -bed'

    _keep_all = ' '
    if keep_all:
        _keep_all = ' -keepAll'

    cmd = f'scanMotifGenomeWide.pl {motif_file_path} {genome}{_bed_out}{_keep_all} -p {cpu} -int'

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


def check_bed_nonoverlap(bed_path):
    cur_chrom = 'NOTACHROM'
    last_end = -1
    line_number = 0
    with open(bed_path) as f:
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


def get_tabix_chroms(file_path):
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
    allc_chroms = get_tabix_chroms(allc_path)
    bed_chroms = get_tabix_chroms(bed_path)
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


def batch_get_motif_rate(cell_allc_table, bed_path, chrom_size_path,
                         out_path, mc_type, cpu=10):
    cell_allc_paths = pd.read_table(cell_allc_table, index_col=0, squeeze=True, header=None)
    if not isinstance(cell_allc_paths, pd.Series):
        raise ValueError('cell_allc_table is malformed. It should only contain 2 columns without header line, '
                         'the 1st is cell id, the 2nd is the allc path.')

    with ProcessPoolExecutor(max_workers=cpu) as executor:
        future_result = {executor.submit(_get_non_overlap_motif_rate,
                                         allc_path=allc_path,
                                         bed_path=bed_path,
                                         chrom_size_path=chrom_size_path,
                                         mc_type=mc_type): cell_id
                         for cell_id, allc_path in cell_allc_paths.iteritems()}
        records = []
        for future in as_completed(future_result):
            cell_id = future_result[future]
            try:
                result = future.result()
                result['cell_id'] = cell_id
                records.append(result)
            except Exception as exc:
                log.info('%r generated an exception: %s' % (cell_id, exc))
    total_result = pd.concat(records, sort=True, ignore_index=True)
    total_result.to_csv(out_path, sep='\t', index=None)
    return
