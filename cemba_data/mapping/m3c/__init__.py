import gzip
import re
import subprocess

import dnaio
import pandas as pd
import pysam


def split_fastq_reads(fastq_path, output_path, trim_b=0, size_l=40, size_r=40, size_m=30):
    """
    Split reads in the fastq file into three parts for remapping.
    Depending on the read length, reads may be
    1) skipped,
    2) split into left and right parts,
    3) split into left, right and middle parts

    left size size_l (name with -l suffix)
    right size size_r (name with -r suffix)
    middle size size_m (name with -m suffix)

    Parameters
    ----------
    fastq_path
    output_path
    trim_b
    size_l
    size_r
    size_m

    Returns
    -------

    """
    trim_b = int(trim_b)
    size_max = max(size_l, size_r)
    with dnaio.open(fastq_path) as f, \
            dnaio.open(output_path, mode='w') as out_f:
        for read in f:
            if trim_b > 0:
                read = read[trim_b:-trim_b]
            read_length = len(read)
            if read_length <= size_max:
                continue
            else:
                # split reads to left and right part, they may have overlap
                left_read = read[:size_l]
                left_read.name += '-l'
                out_f.write(left_read)

                right_read = read[-size_r:]
                right_read.name += '-r'
                out_f.write(right_read)

                # if the middle part is longer enough, we also use it
                if read_length >= (size_l + size_r + size_m):
                    middle_read = read[size_l:-size_r]
                    middle_read.name += '-m'
                    out_f.write(middle_read)
    return


def generate_contacts(bam_path, output_path, chrom_size_path, min_gap=1000):
    """
    Assume the input bam file is sort by name,
    1) generate a split table, contain eight columns, corresponding to split parts of R1 and R2
    2) remove duplicates of the split table
    3) generate contacts based on the split table

    Parameters
    ----------
    bam_path
    output_path
    chrom_size_path
    min_gap

    Returns
    -------

    """
    splits = ['1', '1-l', '1-m', '1-r', '2-l', '2-m', '2-r', '2']
    split_type_pattern = re.compile('-[lrm]$')
    split_table_path = f'{output_path}.split.tsv.unsort'

    # check_seq=False, otherwise pysam raise error for empty bam
    with pysam.AlignmentFile(bam_path, check_sq=False) as bam, \
            open(split_table_path, 'w') as out_f:
        # init
        pre_name_base = ''
        record = {s: '' for s in splits}

        for read in bam:
            name = read.query_name
            name_base = name.split('_')[0]
            # pysam always think they are unpaired cause this is SE bam, so I determine manually
            read_type = name.split('_')[1][0]  # read_type is '1' or '2'
            map_to_reverse_strand = read.flag & 16
            split_type = split_type_pattern.findall(name)
            split_type = split_type[0] if len(split_type) == 1 else ''
            split = f'{read_type}{split_type}'

            # reset if name_base changed
            if name_base != pre_name_base:
                out_f.write('\t'.join([record[s] for s in splits]) + '\n')
                pre_name_base = name_base
                record = {s: '' for s in splits}

            if map_to_reverse_strand:
                if read_type == '1':
                    strand = 0
                else:
                    strand = 1
            else:
                if read_type == '1':
                    strand = 1
                else:
                    strand = 0

            if strand:
                pos = f'{strand}:{read.reference_name}:{read.pos + 1}'
            else:
                pos = f'{strand}:{read.reference_name}:{read.pos + read.rlen}'
            record[split] = pos

        # last record
        out_f.write('\t'.join([record[s] for s in splits]) + '\n')

    # remove duplicates by sorting all fields
    cmd = f'sort -k1 -u {split_table_path} -o {split_table_path[:-7]} && rm -f {split_table_path}'
    subprocess.run(cmd, check=True, shell=True)
    split_table_path = split_table_path[:-7]
    _parse_split_table(split_table_path, output_path, chrom_size_path, min_gap=min_gap)
    return


def _parse_split_table(input_path, output_path, chrom_size_path, min_gap=1000):
    """Generate contacts from the split table"""
    chrom_list = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0).index.tolist()
    chrom_set = set(chrom_list)
    print('Use these chromosomes in contact file:', chrom_set)
    total = 0
    cis_long = 0
    cis_short = 0
    trans = 0

    # chunk reading reduce memory usage
    data = pd.read_csv(input_path, header=0, sep='\t', chunksize=5000)
    with gzip.open(output_path, 'wt') as f:
        for chunk in data:
            for _, row in chunk.iterrows():
                total += 1
                pos_to_chr = {p: p.split(':')[1]
                              for p in row.dropna()
                              if p.split(':')[1] in chrom_set}
                if len(pos_to_chr) < 2:
                    # only one frag mapped, not contact
                    continue
                else:
                    unique_chroms = set(pos_to_chr.values())
                    if len(unique_chroms) > 2:
                        # more than two chr is abnormal
                        continue
                    elif len(unique_chroms) == 2:
                        # trans contact
                        trans += 1
                        # order pos by chrom list, take one from each side/chrom
                        (read1, _), *_, (read2, _) = sorted(pos_to_chr.items(),
                                                            key=lambda i: chrom_list.index(i[1]))
                        f.write(f'{read1}\t0\t{read2}\t1\n'.replace(':', '\t'))
                    else:
                        # cis contact
                        # order the pos, take left or right most
                        read1, *_, read2 = sorted(pos_to_chr.keys(),
                                                  key=lambda i: int(i.split(':')[-1]))
                        pos1 = int(read1.split(':')[-1])
                        pos2 = int(read2.split(':')[-1])
                        if abs(pos1 - pos2) > min_gap:
                            cis_long += 1
                            f.write(f'{read1}\t0\t{read2}\t1\n'.replace(':', '\t'))
                        else:
                            cis_short += 1
    counts = {'CisShortContact': cis_short,
              'CisLongContact': cis_long,
              'TransContact': trans}
    pd.Series(counts).to_csv(f'{output_path}.counts.txt', header=False)
    return
