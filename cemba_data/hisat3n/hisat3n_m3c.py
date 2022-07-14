import re
import shlex
import subprocess
from collections import defaultdict
from itertools import combinations

import dnaio
import pandas as pd
import pysam

# R1 is G to A mutated in unmethylated C
R1_CUT_SITES = {
    'CATG',  # NlaIII
    'CATA',  # NlaIII
    'GATC',  # DpnII or MboI
    'AATC'  # DpnII or MboI
}

# R2 is C to T mutated in unmethylated C
R2_CUT_SITES = {
    'CATG',  # NlaIII
    'TATG',  # NlaIII
    'GATC',  # DpnII or MboI
    'GATT'  # DpnII or MboI
}


def _span_combination(spans, min_length):
    """
    Turn spans into a list of slice object, each slice is a combination of spans

    For example, if spans is [(3, 7), (7, 11), (68, 72), (144, 148)], and min_length is 30,
    the slices will be: [
        slice(3, 72, None),
        slice(3, 148, None),
        slice(7, 72, None),
        slice(7, 148, None),
        slice(68, 148, None)
    ]
    """
    slices = []
    if len(spans) == 0:
        raise
    elif len(spans) == 1:
        base_start, base_end = spans[0]
        slices.append(slice(base_start, base_end))
    else:
        n_span = len(spans)
        for length in range(1, n_span + 1):
            for start in range(0, n_span - length):
                base_start = spans[start][0]  # start of start
                base_end = spans[start + length][1]  # end of end
                if base_end - base_start >= min_length:
                    # remove small slices
                    slices.append(slice(base_start, base_end))
    return slices


def _make_split_pattern(read_type):
    """Get read split pattern based on read type"""
    if str(read_type)[-1] == '1':
        _alt_combine = '|'.join(R1_CUT_SITES)
        # like r'(CATG|GATC|CATA|AATC)'
        split_pattern = re.compile(f'({_alt_combine})')
    elif str(read_type)[-1] == '2':
        _alt_combine = '|'.join(R2_CUT_SITES)
        split_pattern = re.compile(f'({_alt_combine})')
    else:
        raise ValueError(f'read_type must be R1 or R2, got {read_type}')
    return split_pattern


def _split_read_and_make_combination(read, split_pattern, min_length=30):
    """
    Split read by enzyme cut sites and make combination of all possible cut sites

    Parameters
    ----------
    read
        Input read
    split_pattern
        Regex pattern to split read
    min_length
        Minimum read length to keep

    Returns
    -------
    read_split_iter
    """
    # search enzyme cut sites
    read_length = len(read.sequence)
    spans = [i.span() for i in re.finditer(split_pattern, read.sequence)]
    if len(spans) == 0:
        # if there is no cutsite
        spans = [(0, 0), (read_length, read_length)]
    else:
        if spans[0][0] != 0:
            # if the first span does not start at 0, add (0, 0) to capture from the start
            spans = [(0, 0)] + spans
        if spans[-1][1] != read_length:
            # if the last span not end at read_length, add (read_length, read_length) to capture to the end
            spans.append((read_length, read_length))
    # example spans: [(0, 0), (3, 7), (7, 11), (68, 72), (132, 132)]

    # split reads
    slices = _span_combination(spans, min_length=min_length)
    for _slice in slices:
        _read = read[_slice]
        yield _read, _slice


def _trim_site(read_and_slice, read_type):
    """
    Remove DpnII or MboI site from the left read, the site will be included in the right read;
    remove NalIII site from the right read, the site will be included in the left read.
    """
    read, read_slice = read_and_slice
    start = read_slice.start
    stop = read_slice.stop

    sequence = read.sequence
    if read_type == 'R1':
        if sequence[-3:] == 'ATC':
            # If 3' is DpnII or MboI site, clip it (from the left read)
            read = read[:-4]
            stop -= 4
        if sequence[:3] == 'CAT':
            # If 5' is NalIII site, clip it (from the right read)
            read = read[4:]
            start += 4
    else:
        if sequence[-4:-1] == 'GAT':
            # If 3' is DpnII or MboI site, clip it (from the left read)
            read = read[:-4]
            stop -= 4
        if sequence[1:4] == 'ATG':
            # If 5' is NalIII site, clip it (from the right read)
            read = read[4:]
            start += 4
    read.name += f':{start}:{stop}'
    return read, (start, stop)


def split_hisat3n_unmapped_reads(fastq_path,
                                 output_path,
                                 min_length=30):
    """
    Split trimmed fastq file by all possible enzyme cut sites and save to a new fastq file
    Read name is modified, last two field separated by ":" indicate the span start and stop
    A single read might have multiple way of cut site combination, resulting in multiple overlapped reads.
    These reads will all be aligned in single-end mapping, and the post-alignment process
    will determine the unique way of cut site combination.

    Parameters
    ----------
    fastq_path
        Input fastq path
    output_path
        Output fastq path
    min_length
        Minimum read length to keep

    Returns
    -------

    """
    # split pattern
    r1_split_pattern = _make_split_pattern('R1')
    r2_split_pattern = _make_split_pattern('R2')

    with dnaio.open(fastq_path) as f, \
            dnaio.open(output_path, mode='w') as f_out:
        for read in f:
            # read type and split pattern
            read_type = read.name.replace('_', ' ').split(' ')[1][0]
            split_pattern = r1_split_pattern if read_type == '1' else r2_split_pattern

            read_split_iter = _split_read_and_make_combination(
                read=read, split_pattern=split_pattern, min_length=min_length)
            range_set = set()
            for i, read_and_slice in enumerate(read_split_iter):
                # remove overlapping enzyme cut site
                read_split, read_range = _trim_site(read_and_slice, read_type)
                if read_range in range_set:
                    # in some rare cases, the cut site is next to each other,
                    # trim site might cause duplicated ranges and therefore duplicated read names
                    # which is not allowed in the aligners
                    continue
                else:
                    range_set.add(read_range)
                    f_out.write(read_split)
    return


class ReadOverlapGroup:
    """
    Collect overlapped reads within certain span based on the genome coordinates
    This class is used to collect reads aligned in close genome locations
    """

    def __init__(self, read, span=0):
        # genomic position of the read group
        self.chrom = read.reference_name
        self.start = read.reference_start
        self.end = read.reference_end

        # extend region
        self.start -= span
        self.end += span
        self.span = span
        # span apply to both cur and new read, so span 1000 means 2000 distance

        self.reads = [read]

    def add_if_overlap(self, read):
        if read.reference_name != self.chrom:
            # not the same chrom
            return False

        # apply span to new read
        start = read.reference_start
        end = read.reference_end
        start -= self.span
        end += self.span

        if (end > self.start) and (start < self.end):
            self.reads.append(read)
            self.start = min(self.start, start)
            self.end = max(self.end, end)
            # overlap
            return True
        else:
            # not overlap
            return False


class ReadSplitOverlapGroup:
    """
    Collect overlapped reads based on read cut site split position
    This class is used to collect overlapping read parts, genome coordinates are not used
    """

    def __init__(self, read):
        # genomic position of the read group
        self.start = read.get_tag('SS')  # SS is the start position of the read slice
        self.end = read.get_tag('SE')  # SE is the end position of the read slice
        self.reads = [read]
        self.is_read1 = read.is_read1

    def add_if_overlap(self, read):
        if read.is_read1 != self.is_read1:
            # not the same read type
            return False

        # apply span to new read
        start = read.get_tag('SS')
        end = read.get_tag('SE')

        if (end > self.start) and (start < self.end):
            self.reads.append(read)
            self.start = min(self.start, start)
            self.end = max(self.end, end)
            # overlap
            return True
        else:
            # not overlap
            return False


def _remove_overlapped_split_read_parts(read_parts):
    """
    Due to the split read step, some read parts are overlapped.
    This function first calculate the overlapped read part groups,
    then select the longest read part from each group.
    This step does not consider where the reads in each group map to,
    it only considers the relevant relationship within the read coordinates using SS and SE tag.

    read_parts:
    R1 part 1: ----------------
    R1 part 2:        -----------  # this will be removed
    R1 part 3:                         -----------
    R2 part 1: -----------
    R2 part 2:                -----------

    final_parts
    R1 part 1: ----------------
    R1 part 3:                         -----------
    R2 part 1: -----------
    R2 part 2:                -----------
    # Note: R1 and R2 parts are always considered non-overlapping

    Parameters
    ----------
    read_parts
        Read parts from the same pair of reads.

    Returns
    -------
    final non-overlapping read parts
    """
    # determine overlap read groups
    read_groups = []
    for read in read_parts:
        added = False
        for rg in read_groups:
            # added will be True if read overlap with existing rg
            added = rg.add_if_overlap(read)
            if added:
                break
        if not added:
            # no overlap to existing rg, create a new one
            read_groups.append(ReadSplitOverlapGroup(read))

    # select the longest alen read from each group
    final_reads = []
    for group in read_groups:
        longest_read = sorted(
            group.reads,
            key=lambda r: r.get_tag('SE') - r.get_tag('SS'),
            reverse=True)[0]
        final_reads.append(longest_read)
    return final_reads


def _two_read_contact_type(read_1, read_2):
    """
    Determine the contact type between two reads. The contact type includes:
    1. ciscut: the two reads are split from the same read at the cut site, and map to the same chromosome
    2. transcut: the two reads are split from the same read at the cut site, and map to different chromosomes
    3. chimeric: the two reads are split from the same read, but not at the cut site, this might be due to
                 artificial chimeric synthesis event
    4. cis: the two reads are from different read types, and map to the same chromosome
    5. trans: the two reads are from different read types, and map to different chromosomes

    Parameters
    ----------
    read_1
        Read 1
    read_2
        Read 2
    Returns
    -------
    contact type str
    """
    is_same_read = read_1.is_read1 == read_2.is_read1
    # determine if these two reads are split from one original read at the cut site
    if is_same_read:
        is_cutsite = (read_1.get_tag('SS') == read_2.get_tag('SE')) or \
                     (read_1.get_tag('SE') == read_2.get_tag('SS'))
    else:
        is_cutsite = False

    # determine if this is cis contact or trans contact
    cis = read_1.reference_name == read_2.reference_name

    # determine final contact type
    if is_same_read:
        if is_cutsite:
            return 'ciscut' if cis else 'transcut'
        else:
            return 'chimeric'
    else:
        return 'cis' if cis else 'trans'


def _extract_contact_info(reads, span=1000):
    """
    Extract chromatin contacts from a list of aligned reads

    Parameters
    ----------
    reads
        aligned reads
    span
        reads within this distance will be treated as single fragment and merged before consider contacts.

    Returns
    -------
    list of possible chromatin contacts with contact type judgement
    """
    # merge reads with a span, ignore strand or read type difference
    read_groups = []
    for read in reads:
        added = False
        for rg in read_groups:
            # added will be True if read overlap with existing rg
            added = rg.add_if_overlap(read)
            if added:
                break
        if not added:
            # no overlap to existing rg, create a new one
            read_groups.append(
                ReadOverlapGroup(read,
                                 span=span))

    n_groups = len(read_groups)
    if n_groups < 2:
        # not a contact (after consider span)
        return [(None, 'no')]
    elif n_groups >= 2:
        results = []
        multi = '_multi' if n_groups > 2 else ''  # one read pair contains multiple read contacts
        for g1, g2 in combinations(read_groups, 2):
            group_contacts = {}
            # read_1 and read_2 are the two reads in different groups, not necessarily correspond to R1 or R2
            for read_1 in g1.reads:
                for read_2 in g2.reads:
                    contact_type = _two_read_contact_type(read_1, read_2)
                    group_contacts[contact_type] = (read_1, read_2)
            # if g1 and g2 both containing multiple reads, we will only choose one contact type
            # to represent the relationship between the two groups
            # the choice is ordered by the following priority:
            if 'ciscut' in group_contacts:
                results.append((group_contacts['ciscut'], 'ciscut' + multi))
            elif 'transcut' in group_contacts:
                results.append((group_contacts['transcut'], 'transcut' + multi))
            elif 'chimeric' in group_contacts:
                # if chimeric, we will not consider the contact,
                # saving the contacts here is for quantification purpose later.
                results.append((None, 'chimeric'))
            elif 'cis' in group_contacts:
                results.append((group_contacts['cis'], 'cis' + multi))
            elif 'trans' in group_contacts:
                results.append((group_contacts['trans'], 'trans' + multi))
        return results


class ContactWriter:
    def __init__(self, output_prefix):
        self.output_path = f'{output_prefix}.raw_contacts.tsv'
        self.counter = defaultdict(int)
        self.cur_read_pair = 0

    def __enter__(self):
        self.out = open(self.output_path, 'w')
        title = 'read_pair_id\tcontact_type\t' \
                'chrom1\tstart1\tend1\tstrand1\tread_start1\tread_end1\tread_type1\t' \
                'chrom2\tstart2\tend2\tstrand2\tread_start2\tread_end2\tread_type2\n'
        self.out.write(title)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.out.close()

    @classmethod
    def _read_to_string_contact(cls, read):
        read_type = 1 if read.is_read1 else 2
        return f'{read.reference_name}\t{read.pos}\t{read.aend}\t' \
               f'{int(read.is_reverse)}\t{read.get_tag("SS")}\t{read.get_tag("SE")}\t{read_type}'

    def write(self, results):
        for result in results:
            read_pair, contact_type = result
            self.counter[contact_type] += 1
            if read_pair is None:
                # this happens when contact type is chimeric, we will not save the contact
                continue
            read_1, read_2 = read_pair
            read_1_str = self._read_to_string_contact(read_1)
            read_2_str = self._read_to_string_contact(read_2)
            contact_row = f"{self.cur_read_pair}\t{contact_type}\t{read_1_str}\t{read_2_str}\n"
            self.out.write(contact_row)
        self.cur_read_pair += 1
        return


def _dedup_chrom_df(chrom_df):
    """Deduplicate contacts within the same chromosome"""
    # determine duplicated coords for each column
    dup_judge = chrom_df[['start1', 'start2', 'end1', 'end2']].apply(lambda i: i.duplicated())

    # if contact on both sides have one position (start or end) being the same,
    # it is considered as a duplicated contact
    dup_call = (dup_judge['start1'] | dup_judge['end1']) & (dup_judge['start2'] | dup_judge['end2'])
    chrom_df = chrom_df[~dup_call]
    return chrom_df


def _dedup_contacts(output_prefix):
    """Deduplicate contacts by chromosome and position"""
    input_path = f'{output_prefix}.raw_contacts.tsv'
    output_path = f'{output_prefix}.dedup_contacts.tsv.gz'

    contacts = pd.read_csv(input_path, sep='\t')
    contacts = contacts.sort_values(
        by=['chrom1', 'chrom2', 'start1', 'start2', 'end1', 'end2'])

    total_dedup = []
    for _, chrom_df in contacts.groupby(['chrom1', 'chrom2']):
        total_dedup.append(_dedup_chrom_df(chrom_df))
    total_dedup = pd.concat(total_dedup)
    total_dedup.to_csv(output_path, sep='\t', index=False)
    subprocess.run(shlex.split(f'rm -f {input_path}'), check=True)
    return


def _contact_to_hic_format(output_prefix):
    contact_df = pd.read_csv(f'{output_prefix}.dedup_contacts.tsv.gz', sep='\t')
    # some columns used in hic format
    contact_df['hic0'] = 0
    contact_df['hic1'] = 1
    output_path = f'{output_prefix}.3C.contact.tsv.gz'
    contact_df[[
        'strand1', 'chrom1', 'start1', 'hic0', 'strand2', 'chrom2', 'start2',
        'hic1'
    ]].to_csv(output_path, header=None, index=None, sep='\t')
    return


def remove_overlap_read_parts(in_bam_path, out_bam_path):
    with pysam.AlignmentFile(in_bam_path) as bam, \
            pysam.AlignmentFile(out_bam_path, header=bam.header, mode='w') as out_bam:
        count = 0
        cur_read_pair_name = None
        cur_read_parts = []
        for read in bam:
            read_pair_name, others = read.qname.split('_')
            read_type, start, stop = others.split(':')
            read.is_read1 = read_type == '1'
            read.is_read2 = not read.is_read1
            read.set_tag('SS', int(start))  # SS is the start position of read slice
            read.set_tag('SE', int(stop))  # SE is the end position of read slice
            if read_pair_name == cur_read_pair_name:
                cur_read_parts.append(read)
            else:
                # read to a new read pair
                if len(cur_read_parts) > 0:
                    # process the previous read pair
                    count += 1
                    final_reads = _remove_overlapped_split_read_parts(cur_read_parts)
                    for final_read in final_reads:
                        # need to put back the flag
                        # otherwise picard will raise errors
                        # dedup still perform at single read level
                        # unless the flag can be properly set
                        final_read.is_read1 = False
                        final_read.is_read2 = False
                        out_bam.write(final_read)
                # initiate the next pair
                cur_read_pair_name = read_pair_name
                cur_read_parts = [read]
        if len(cur_read_parts) > 0:
            # process the last read pair
            final_reads = _remove_overlapped_split_read_parts(cur_read_parts)
            for final_read in final_reads:
                final_read.is_read1 = False
                final_read.is_read2 = False
                out_bam.write(final_read)
    return


def call_chromatin_contacts(bam_path: str,
                            contact_prefix: str,
                            dedup_contact: bool = True,
                            save_hic_format: bool = True):
    """
    Process 3C bam file and generate contact file.

    Parameters
    ----------
    bam_path : str
        Path to 3C bam file.
    contact_prefix: str
        Prefix of output contact file.
    dedup_contact : bool, optional
        Whether to deduplicate contacts.
    save_hic_format : bool, optional
        Whether to save the contact in hic format.

    Returns
    -------

    """
    with pysam.AlignmentFile(bam_path) as bam, \
            ContactWriter(contact_prefix) as out_contacts:
        count = 0
        cur_read_pair_name = None
        cur_read_parts = []
        for read in bam:
            read_pair_name = read.qname.split('_')[0]
            if not read.has_tag('SS'):
                read.set_tag('SS', 0)
                read.set_tag('SE', read.qlen)

            if read_pair_name == cur_read_pair_name:
                cur_read_parts.append(read)
            else:
                # read to a new read pair
                if len(cur_read_parts) > 0:
                    # process the previous read pair
                    count += 1
                    results = _extract_contact_info(cur_read_parts, span=1000)
                    out_contacts.write(results)
                # initiate the next pair
                cur_read_pair_name = read_pair_name
                cur_read_parts = [read]
        if len(cur_read_parts) > 0:
            # process the last read pair
            results = _extract_contact_info(cur_read_parts, span=1000)
            out_contacts.write(results)

    if dedup_contact:
        # dedup contacts
        _dedup_contacts(contact_prefix)

    if save_hic_format:
        # save hic format
        _contact_to_hic_format(contact_prefix)

    # save counts
    stat = out_contacts.counter
    stat['mapped_frag'] = count
    pd.Series(stat).to_csv(f'{contact_prefix}.contact_stats.csv', header=False)
    return
