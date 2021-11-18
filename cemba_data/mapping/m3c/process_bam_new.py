import pysam
import subprocess
from itertools import combinations
from collections import defaultdict
import pandas as pd


class ReadOverlapGroup:
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
    def __init__(self, read):
        # genomic position of the read group
        self.start = read.get_tag('SS')
        self.end = read.get_tag('SE')
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


def process_bam_read_pairs(read_parts):
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


def two_read_contact_type(read_1, read_2):
    is_same_read = read_1.is_read1 == read_2.is_read1
    if is_same_read:
        is_cutsite = (read_1.get_tag('SS') == read_2.get_tag('SE')) or \
                     (read_1.get_tag('SE') == read_2.get_tag('SS'))
    else:
        is_cutsite = False

    cis = read_1.reference_name == read_2.reference_name

    if is_same_read:
        if is_cutsite:
            return 'ciscut' if cis else 'transcut'
        else:
            return 'chimeric'
    else:
        return 'cis' if cis else 'trans'


def extract_contact_info(reads, span=1000):
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
        multi = '_multi' if n_groups > 2 else ''
        for g1, g2 in combinations(read_groups, 2):
            group_contacts = {}
            for read_1 in g1.reads:
                for read_2 in g2.reads:
                    contact_type = two_read_contact_type(read_1, read_2)
                    group_contacts[contact_type] = (read_1, read_2)
            if 'ciscut' in group_contacts:
                results += [(group_contacts['ciscut'], 'ciscut' + multi)]
            elif 'transcut' in group_contacts:
                results += [(group_contacts['transcut'], 'transcut' + multi)]
            elif 'chimeric' in group_contacts:
                results += [(None, 'chimeric')]
            elif 'cis' in group_contacts:
                results += [(group_contacts['cis'], 'cis' + multi)]
            elif 'trans' in group_contacts:
                results += [(group_contacts['trans'], 'trans' + multi)]
        return results


def _read_to_string_contact(read):
    read_type = 1 if read.is_read1 else 2
    return f'{read.reference_name}\t{read.pos}\t{read.aend}\t' \
           f'{int(read.is_reverse)}\t{read.get_tag("SS")}\t{read.get_tag("SE")}\t{read_type}'


class ContactWriter:
    def __init__(self, output_prefix):
        self.out = open(f'{output_prefix}.raw_contacts.tsv', 'w')
        title = 'read_pair_id\tcontact_type\t' \
                'chrom1\tstart1\tend1\tstrand1\tread_start1\tread_end1\tread_type1\t' \
                'chrom2\tstart2\tend2\tstrand2\tread_start2\tread_end2\tread_type2\n'
        self.out.write(title)
        self.counter = defaultdict(int)
        self.cur_read_pair = 0

    def dump(self, results):
        for result in results:
            read_pair, contact_type = result
            self.counter[contact_type] += 1
            if read_pair is None:
                continue
            read_1, read_2 = read_pair
            read_1_str = _read_to_string_contact(read_1)
            read_2_str = _read_to_string_contact(read_2)
            contact_row = f"{self.cur_read_pair}\t{contact_type}\t{read_1_str}\t{read_2_str}\n"
            self.out.write(contact_row)
        self.cur_read_pair += 1
        return

    def close(self):
        self.out.close()


def _dedup_chrom_df(chrom_df):
    dup_judge = chrom_df[['start1', 'start2', 'end1', 'end2']].apply(lambda i: i.duplicated())
    # if reads on both sides have one postion (start or end) being the same,
    # it is consider as a duplciated contact
    dup_call = (dup_judge['start1'] | dup_judge['end1']) & (dup_judge['start2'] | dup_judge['end2'])
    chrom_df = chrom_df[~dup_call]
    return chrom_df


def _dedup_contacts(output_prefix):
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
    subprocess.run(f'rm -f {input_path}', shell=True)
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


def process_3c_bam(bam_path, bam_prefix, contact_prefix):
    bam = pysam.AlignmentFile(bam_path)

    out_bam = pysam.AlignmentFile(f'{bam_prefix}.non_overlap_read_parts.bam',
                                  header=bam.header,
                                  mode='w')
    out_contacts = ContactWriter(contact_prefix)
    count = 0
    cur_read_pair_name = None
    cur_read_parts = []
    for read in bam:
        read_pair_name, others = read.qname.split('_')
        read_type, *_, start, stop = others.split(':')
        read.is_read1 = read_type == '1'
        read.is_read2 = not read.is_read1
        read.set_tag('SS', int(start))
        read.set_tag('SE', int(stop))
        if read_pair_name == cur_read_pair_name:
            cur_read_parts.append(read)
        else:
            # read to a new read pair
            if len(cur_read_parts) > 0:
                # process the previous read pair
                count += 1
                final_reads = process_bam_read_pairs(cur_read_parts)
                results = extract_contact_info(final_reads, span=1000)
                out_contacts.dump(results)
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
        final_reads = process_bam_read_pairs(cur_read_parts)
        results = extract_contact_info(final_reads, span=1000)
        out_contacts.dump(results)
        for final_read in final_reads:
            # need to put back the flag
            # otherwise picard will raise errors
            # dedup still perform at single read level
            # unless the flag can be properly set
            final_read.is_read1 = False
            final_read.is_read2 = False
            out_bam.write(final_read)
    out_bam.close()
    out_contacts.close()

    # dedup contacts
    _dedup_contacts(contact_prefix)
    # save hic format
    _contact_to_hic_format(contact_prefix)

    # save counts
    stat = out_contacts.counter
    stat['mapped_frag'] = count
    pd.Series(stat).to_csv(f'{contact_prefix}.3c_bam_stats.csv', header=False)
    return
