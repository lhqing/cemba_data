from heapq import heappop, heappush
import numpy as np
from cemba_data.tools.allc import extract_allc
import subprocess
from ._open import open_gz, open_allc
from ..utilities import parse_mc_pattern, parse_chrom_size


class _Region:
    def __init__(self, line):
        try:
            self.chrom, start, end, self.uid = line.strip().split('\t')
            self.start = int(start)
            self.end = int(end)
        except ValueError as e:
            if line == '':
                # end of file
                raise EOFError
            else:
                # file malformed
                raise e
        return

    def __lt__(self, region):
        """implement < for heappushpop to compare"""
        return self.end < region.end

    def __gt__(self, region):
        """implement > for heappush to compare"""
        return self.end > region.end


class _Site:
    def __init__(self, line):
        try:
            self.chrom, pos, self.strand, self.context, mc, cov, _ = line.strip().split('\t')
            self.pos = int(pos)
            self.mc = int(mc)
            self.cov = int(cov)
        except ValueError as e:
            if line == '':
                # end of file
                raise EOFError
            else:
                # file malformed
                raise e
        return


class _CurrentQueue:
    """Keep a heap order by region.end"""

    def __init__(self):
        self.heap = []
        self.record_dict = {}

        # use temp to prevent for loop for every site
        # only restore temp to record_dict when reach to region end
        self.temp_mc = 0
        self.temp_cov = 0
        return

    def add_region(self, region):
        self.record_dict[region.uid] = np.zeros(shape=(2,))
        heappush(self.heap, region)

    def dump_temp(self):
        # record temp first
        for value in self.record_dict.values():
            # temp do not contain current site
            value += (self.temp_mc, self.temp_cov)
        self.temp_cov = 0
        self.temp_mc = 0

    def check_site(self, site, new_regions):
        return_content = []

        # pop region until the end > pos
        while (len(self.heap) > 0) and (site.pos > self.heap[0].end):
            if self.temp_cov != 0:
                self.dump_temp()
            _current_region = heappop(self.heap)
            # deal with region end, use while because multiple region could end in same pos
            return_count = self.record_dict.pop(_current_region.uid).astype(int)

            return_line = [_current_region.chrom, _current_region.start,
                           _current_region.end, _current_region.uid,
                           return_count[0], return_count[1]]
            return_content.append(return_line)

        # add new region next
        for new_region in new_regions:
            # check already finished region before add
            if site.pos > new_region.end:
                return_line = [new_region.chrom, new_region.start,
                               new_region.end, new_region.uid, '.', '.']
                return_content.append(return_line)
            else:
                if self.temp_cov != 0:
                    self.dump_temp()
                self.add_region(new_region)

        # record current site to temp
        self.temp_mc += site.mc  # reset temp and add current site
        self.temp_cov += site.cov  # reset temp and add current site
        return return_content


class _RegionCounter:
    """Only take care single chromosome, single context, and single bed_handle"""

    def __init__(self, bed_handle):
        self.handle = bed_handle
        self.current_queue = _CurrentQueue()

        # init current_queue by 1st region:
        self.next_region = None
        self.next_start = None
        self.update_next()

        self.eof = False
        return

    def update_next(self):
        try:
            self.next_region = _Region(self.handle.readline())
            self.next_start = self.next_region.start
        except EOFError:
            self.eof = True
        return

    def check_site(self, site):
        new_regions = []

        # bed format half open, so > not >=
        while site.pos > self.next_start and not self.eof:
            # use while because there could be multiple region start at same pos
            # put this region into current queue
            new_regions.append(self.next_region)
            self.update_next()

        region_line_list = self.current_queue.check_site(site, new_regions)
        return region_line_list


def _map_to_region_worker(allc_path, bed_path, context_map, mc_patterns, chrom, max_cov_cutoff=2):
    """One bed file, one chromosome, but multiple pattern"""
    allc_handle = open_allc(allc_path, mode='r', region=chrom)
    counter_dict = {}
    record_dict = {}
    for mc_pattern in mc_patterns:
        bed_handle = open_gz(bed_path, mode='r', region=chrom)
        region_counter = _RegionCounter(bed_handle)
        counter_dict[mc_pattern] = region_counter
        record_dict[mc_pattern] = []

    for line in allc_handle:
        try:
            site = _Site(line)
        except EOFError:
            # read to allc end
            break
        if site.cov > max_cov_cutoff:
            continue

        try:
            pattern = context_map[site.context]
            region_counter = counter_dict[pattern]

            records = record_dict[pattern]
            region_line_list = region_counter.check_site(site)
            records.extend(region_line_list)

            if region_counter.eof and len(region_counter.current_queue.record_dict) == 0:
                counter = counter_dict.pop(pattern)
                counter.handle.close()

        except KeyError:
            # context not included in pattern
            continue

    # read the remaining row in bed, count should all be 0 because allc already finished
    for mc_pattern, counter in counter_dict.items():
        records = record_dict[mc_pattern]
        fake_site = _Site('chrNNN\t9999999999\t/\t---\t0\t0\t1')
        # get all remaining regions
        region_line_list = counter.check_site(fake_site)
        records.extend(region_line_list)

    # close handles
    allc_handle.close()
    for counter in counter_dict.values():
        counter.handle.close()

    for records in record_dict.values():
        records.sort(key=lambda i: i[1])  # inplace sort by region start

    return record_dict


def _map_to_region_non_overlap_worker(allc_path, chrom, bed_path, context_map, mc_patterns):
    """One non-overlap bed file, one chromosome, but multiple pattern"""
    record_dict = {pattern: [] for pattern in mc_patterns}
    temp_dict = {pattern: np.zeros(shape=(2,), dtype=np.uint32)
                 for pattern in mc_patterns}

    with open_allc(allc_path, mode='r', region=chrom) as f, open_gz(bed_path, mode='r', region=chrom) as bed:
        # TODO empty bed
        cur_region = _Region(bed.readline())
        for line in f:
            site = _Site(line)
            pattern = context_map[site.context]
            if site.pos <= cur_region.start:
                continue
            elif site.pos <= cur_region.end:
                # add region
                temp_dict[pattern] += (site.mc, site.cov)
            else:
                # restore region and reset
                for pattern in mc_patterns:
                    region_line = [cur_region.chrom, cur_region.start, cur_region.end, cur_region.uid,
                                   temp_dict[pattern][0], temp_dict[pattern][1]]
                    record_dict[pattern].append(region_line)
                cur_region = _Region(bed.readline())
                # TODO zero region

                temp_dict = {pattern: np.array(shape=(2,), dtype=np.uint32)
                             for pattern in mc_patterns}
    return


def map_to_region(allc_path, bed_paths, bed_names, mc_patterns,
                  chromsize_path, out_prefix, max_cov_cutoff=2):
    context_map = {}
    for mc_pattern in mc_patterns:
        context_set = parse_mc_pattern(mc_pattern)
        for context in context_set:
            if context in context_map:
                raise ValueError('mc_patterns contain mC patterns that overlap. '
                                 'Each map-to-region only take exclusive mC pattern. '
                                 'Run map-to-region multiple times for this situation.')
            context_map[context] = mc_pattern
    out_handle_dict = {}
    for mc_pattern in mc_patterns:
        for bed_name in bed_names:
            out_path = out_prefix + f'{bed_name}_{mc_pattern}.count_table.bed.gz'
            out_handle = open_gz(out_path, mode='w')
            out_handle_dict[(mc_pattern, bed_name)] = out_handle

    chrom_dict = parse_chrom_size(chromsize_path, remove_chr_list=None)

    # TODO: this can be paralleled
    for chrom in chrom_dict:
        for bed_name, bed_path in zip(bed_names, bed_paths):
            record_dict = _map_to_region_worker(allc_path, bed_path, context_map, mc_patterns,
                                                chrom, max_cov_cutoff=max_cov_cutoff)
            for mc_pattern, records in record_dict.items():
                out_handle = out_handle_dict[(mc_pattern, bed_name)]
                for line in records:
                    out_handle.write('\t'.join(map(str, line)) + '\n')
    for handle in out_handle_dict.values():
        handle.close()
    return


def _transfer_bin_size(bin_size: int) -> str:
    """Get proper str for a large bin_size"""
    if bin_size > 1000000:
        bin_size_mode = bin_size % 1000000
        bin_size_mode = f'{bin_size_mode/1000000:.1f}'[1:] if bin_size_mode >= 100000 else ''
        bin_size_str = f'{bin_size//1000000}{bin_size_mode}Mb'
    elif bin_size > 1000:
        bin_size_mode = bin_size % 1000
        bin_size_mode = f'{bin_size_mode/1000:.1f}'[1:] if bin_size_mode >= 100 else ''
        bin_size_str = f'{bin_size//1000}{bin_size_mode}Kb'
    else:
        bin_size_str = f'{bin_size}bp'
    return bin_size_str


"""
Sparse Table Format
#chr1   length_of_chr1
#chr2   length_of_chr2
#chr3   length_of_chr3
#chr4   length_of_chr4
bin_id  bin_mc  bin_cov
...
Used to save small chromosome bins with very small bin size
"""


# TODO: Compare this to the OOP function
def _map_to_sparse_chrom_bin(allc_path, out_prefix, chrom_size_file,
                             remove_additional_chrom=False, bin_size=50):
    """
    Parameters
    ----------
    allc_path
    out_prefix
    chrom_size_file
    remove_additional_chrom
    bin_size

    Returns
    -------

    """
    chrom_dict = parse_chrom_size(chrom_size_file)
    cur_chrom = 'TOTALLY_NOT_A_CHROM'
    cur_chrom_end = 0
    bin_end = min(cur_chrom_end, bin_size)
    bin_id = -1
    temp_mc, temp_cov = 0, 0

    out_prefix = out_prefix.rstrip('.')
    out_path = out_prefix + f'.{_transfer_bin_size(bin_size)}.sparse_table.tsv.gz'

    with open_allc(allc_path) as allc, \
            open_gz(out_path, 'w') as out_handle:
        # add header to indicate chromosome order
        out_handle.write(f'#bin_size\t{bin_size}')
        for chrom, chrom_size in chrom_dict.items():
            out_handle.write(f'#{chrom}\t{chrom_size}\n')

        for line in allc:
            chrom, pos, _, _, mc, cov, *_ = line.split("\t")
            pos = int(pos) - 1  # bed format
            mc = int(mc)
            cov = int(cov)
            if pos >= bin_end or cur_chrom != chrom:
                # write line
                if temp_cov > 0:
                    out_handle.write("\t".join(map(str, [bin_id, temp_mc, temp_cov])) + "\n")

                # reset_chrom
                if cur_chrom != chrom:
                    cur_chrom = chrom
                    try:
                        cur_chrom_end = chrom_dict[chrom]
                    except KeyError as e:
                        # chrom not in chrom size file
                        if remove_additional_chrom:
                            continue
                        else:
                            raise e

                # reset bin
                temp_mc, temp_cov = mc, cov
                bin_end = min(cur_chrom_end, pos // bin_size * bin_size + bin_size)
                bin_id += 1
            else:
                temp_mc += mc
                temp_cov += cov
        # write last piece
        if temp_cov > 0:
            out_handle.write("\t".join(map(str, [bin_id, temp_mc, temp_cov])) + "\n")
    return out_path


def map_to_sparse_chrom_bin(allc_path, out_prefix, chrom_size_file,
                            bin_size=500, mc_contexts=None,
                            remove_additional_chrom=False):
    if mc_contexts is not None:
        extract_allc(allc_path=allc_path,
                     out_prefix=out_prefix,
                     mc_contexts=mc_contexts)
        allc_list = []
        for mc_context in mc_contexts:
            extracted_path = out_prefix.rstrip('.') + f'.extract_{mc_context}.tsv.gz'
            allc_list.append(extracted_path)
    else:
        allc_list = [allc_path]

    for _allc_path in allc_list:
        _map_to_sparse_chrom_bin(allc_path=_allc_path,
                                 out_prefix=_allc_path.rstrip('.tsv.gz'),
                                 chrom_size_file=chrom_size_file,
                                 remove_additional_chrom=remove_additional_chrom,
                                 bin_size=bin_size)
        if mc_contexts is not None:
            # this will only remove extracted ALLC, but not the input ALLC
            subprocess.run(['rm', 'f', _allc_path, _allc_path + '.tbi'])
    return


def aggregate_sparse_tables():
    genome_bin_bed = None
    sparse_matrix = None
    pass
    return genome_bin_bed, sparse_matrix
