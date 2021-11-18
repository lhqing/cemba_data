import dnaio
import re

R1_CUT_SITES = {
    'CATG',  # NlaIII
    'CATA',  # NlaIII
    'GATC',  # DpnII or MboI
    'AATC'  # DpnII or MboI
}
R2_CUT_SITES = {
    'CATG',  # NlaIII
    'TATG',  # NlaIII
    'GATC',  # DpnII or MboI
    'GATT'  # DpnII or MboI
}


def _span_combination(spans, min_length):
    """Turn spans into a list of slice object,
    spans is a list such as [(3, 7), (7, 11), (68, 72), (144, 148)]"""
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
    if read_type.upper() == 'R1':
        _alt_combine = '|'.join(R1_CUT_SITES)
        # like r'(CATG|GATC|CATA|AATC)'
        split_pattern = re.compile(f'({_alt_combine})')
    elif read_type.upper() == 'R2':
        _alt_combine = '|'.join(R2_CUT_SITES)
        split_pattern = re.compile(f'({_alt_combine})')
    else:
        raise ValueError(f'read_type must be R1 or R2, got {read_type}')
    return split_pattern


def _split_read_and_make_combination(read, split_pattern, min_length=30):
    # search enzyme cut sites
    read_length = len(read.sequence)
    spans = [i.span() for i in re.finditer(split_pattern, read.sequence)]
    if len(spans) == 0:
        # if there is not cutsite
        spans = [(0, 0), (read_length, read_length)]
    else:
        if spans[0][0] != 0:
            # if the first span not start at 0, add (0, 0) to capture from the start
            spans = [(0, 0)] + spans
        if spans[-1][1] != read_length:
            # if the last span not end at read_length, add (read_length, read_length) to capture to the end
            spans.append([read_length, read_length])
    # example spans: [(0, 0), (3, 7), (7, 11), (68, 72), (132, 132)]

    # split reads
    slices = _span_combination(spans, min_length=min_length)
    for _slice in slices:
        _read = read[_slice]
        yield _read, _slice


def _trim_site(read_and_slice, read_type):
    """Remove DpnII or MboI site from the left read,
    remove NalIII site from the right read"""
    read, slice = read_and_slice
    start = slice.start
    stop = slice.stop

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


def split_reads(fastq_path, output_path, read_type, min_length=30):
    """
    Split trimmed fastq file by all possible enzyme cut sites and save to a new fastq file
    Read name is modified, last two field separated by ":" indicate the span start and stop

    Parameters
    ----------
    fastq_path
        Input fastq path
    output_path
        Output fastq path
    read_type
        R1 or R2
    min_length
        Minimum read length to keep

    Returns
    -------

    """
    # determine split pattern
    split_pattern = _make_split_pattern(read_type)
    with dnaio.open(fastq_path) as f, dnaio.open(output_path,
                                                 mode='w') as out_f:
        for read in f:
            read_split_iter = _split_read_and_make_combination(
                read=read, split_pattern=split_pattern, min_length=min_length)
            range_set = set()
            for i, read_and_slice in enumerate(read_split_iter):
                # remove overlapping enzyme cut site
                read_split, range = _trim_site(read_and_slice, read_type)
                if range in range_set:
                    # in some rare cases, the cut site is next to each other,
                    # trim site might cause duplicated ranges,
                    # which is not allowed in aligners
                    continue
                else:
                    range_set.add(range)
                    out_f.write(read_split)
    return
