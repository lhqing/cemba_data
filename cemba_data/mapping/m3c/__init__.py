import pysam
import subprocess
import pandas as pd
import dnaio


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


def _output(rfh, tot1, tot2, pre_id, locs):
    c1, c2 = locs[:4].count(''), locs[4:].count('')
    if c1 < 4:
        tot1 += 1
    if c2 < 4:
        tot1 += 1
    if (c1 + c2) < 7:
        tot2 += 1
        rfh.write(pre_id + '\t' + '\t'.join(locs) + '\n')
    return tot1, tot2


def _parse_bam(bam_path, output_path):
    splits = ['1', '1-1', '1-3', '1-2', '2-2', '2-3', '2-1', '2']
    split_dict = {x: i for i, x in enumerate(splits)}
    pre_id = ''
    locs = ['' for _ in range(len(splits))]

    dfh = pysam.AlignmentFile(bam_path, 'rb')
    rfh = open(output_path, 'w')
    tot1, tot2 = 0, 0
    for read in dfh:
        line = str(read).split()
        if '_' in line[0]:
            _id = line[0].split('_')[0]
            split_st = line[0].split('_')[1].split(':')[0]
            if line[0][-2:] == '-l':
                split_st += '-1'
            elif line[0][-2:] == '-r':
                split_st += '-2'
            elif line[0][-2:] == '-m':
                split_st += '-3'
            if read.flag & 16:
                if split_st.split('-')[0] == '1':
                    strand = 0
                else:
                    strand = 1
            else:
                if split_st.split('-')[0] == '1':
                    strand = 1
                else:
                    strand = 0
            if _id != pre_id:
                tot1, tot2 = _output(rfh, tot1, tot2, pre_id, locs)
                pre_id = _id
                locs = ['' for _ in range(len(splits))]
            if strand == 1:
                locs[split_dict[split_st]] = f'1:' \
                                             f'{dfh.get_reference_name(read.reference_id)}:' \
                                             f'{str(read.pos + 1)}'
            if strand == 0:
                locs[split_dict[split_st]] = f'0:' \
                                             f'{dfh.get_reference_name(read.reference_id)}:' \
                                             f'{str(read.pos + len(line[9]))}'
    _output(rfh, tot1, tot2, pre_id, locs)
    rfh.close()

    cmd = f'sort -k2 -u {output_path} -o {output_path}'
    subprocess.run(cmd, shell=True, check=True)
    return output_path


def _parse_split_table(input_path, output_path, chrom_size_path, min_gap=2500):
    output_path = str(output_path)
    if output_path.endswith('.gz'):
        # remove gz first
        output_path = output_path[:-3]

    dfh = open(input_path, 'r')
    rfh = open(output_path, 'w')
    chrom = pd.read_csv(chrom_size_path, sep='\t', header=None, index_col=0).index.tolist()
    chr_dict = {x: i for i, x in enumerate(chrom)}
    tot, cis_long, cis_short, trans = 0, 0, 0, 0

    for line in dfh:
        tot += 1
        frag = line.strip('\n').split('\t')[1:]
        i = 0
        read1 = ''
        for i in range(8):
            if frag[i] != '':
                read1 = frag[i].replace(':', '\t')
                break
        for j in range(7, i, -1):
            if frag[j] != '':
                read2 = frag[j].replace(':', '\t')
                read1_chro, read1_pos = read1.split('\t')[1:3]
                read2_chro, read2_pos = read2.split('\t')[1:3]
                read1_pos, read2_pos = int(read1_pos), int(read2_pos)
                if (read1_chro in chrom) and (read2_chro in chrom):
                    if chr_dict[read1_chro] > chr_dict[read2_chro] or (
                            chr_dict[read1_chro] == chr_dict[read2_chro] and read1_pos > read2_pos):
                        read1, read2 = read2, read1
                    if chr_dict[read1_chro] == chr_dict[read2_chro]:
                        if abs(read1_pos - read2_pos) > min_gap:
                            cis_long += 1
                            rfh.write(read1 + '\t0\t' + read2 + '\t1\n')
                        else:
                            cis_short += 1
                    else:
                        trans += 1
                        rfh.write(read1 + '\t0\t' + read2 + '\t1\n')
                break
    dfh.close()
    rfh.close()
    counts = {'CisShortContact': cis_short,
              'CisLongContact': cis_long,
              'TransContact': trans}
    pd.Series(counts).to_csv(f'{output_path}.counts.txt', header=False)

    # sort the contacts
    subprocess.run(f'sort -k2,2 -k6,6 -k1,1n -k5,5n -k3,3n {output_path} -o {output_path}', shell=True, check=True)

    # compress the contacts, .gz will be added
    subprocess.run(f'gzip {output_path}', shell=True, check=True)
    return output_path


def generate_contacts(bam_path, output_path, chrom_size_path, min_gap=2500, keep_split_table=False):
    split_table_path = f'{output_path}.split.tsv'

    # from bam to eight column split table
    split_table_path = _parse_bam(bam_path=bam_path, output_path=split_table_path)

    # from split table to contacts
    _parse_split_table(input_path=split_table_path,
                       output_path=output_path,
                       chrom_size_path=chrom_size_path,
                       min_gap=min_gap)

    if not keep_split_table:
        subprocess.run(['rm', '-f', split_table_path])
    return
