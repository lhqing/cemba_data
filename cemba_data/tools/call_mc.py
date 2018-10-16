import subprocess
import shlex
import pathlib
import gzip
import pandas as pd


def read_faidx(faidx_path):
    return pd.read_table(faidx_path, index_col=0, header=None,
                         names=['NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH'])


def get_chromosome_sequence(fasta_path, fai_df, query_chrom):
    if query_chrom not in fai_df.index:
        return None
    else:
        chrom_pointer = fai_df.loc[query_chrom, 'OFFSET']
        tail = fai_df.loc[query_chrom, 'LINEBASES'] - fai_df.loc[query_chrom, 'LINEWIDTH']
        seq = ""
        with open(fasta_path, 'r') as f:
            f.seek(chrom_pointer)
            for line in f:
                if line[0] == '>':
                    break
                seq += line[:tail]  # trim \n
        return seq


def call_methylated_sites(inputf, reference_fasta,
                          num_upstr_bases=0,
                          num_downstr_bases=2,
                          buffer_line_number=100000,
                          min_mapq=10,
                          path_to_samtools="",
                          min_base_quality=1):
    # Check fasta index
    if not pathlib.Path(reference_fasta + ".fai").exists():
        raise FileNotFoundError("Reference fasta not indexed. Use samtools faidx to index it and run again.")
    fai_df = read_faidx(pathlib.Path(reference_fasta + ".fai"))

    if not pathlib.Path(inputf + ".bai").exists():
        print("Input not indexed. Indexing...")
        subprocess.check_call(shlex.split(path_to_samtools + "samtools index " + inputf))

    # mpileup
    mpileup_cmd = f"samtools mpileup -Q {min_base_quality} -q {min_mapq} -B -f {reference_fasta} {inputf}'"
    pipes = subprocess.Popen(shlex.split(mpileup_cmd),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    result_handle = pipes.stdout

    # Output handel
    input_path = pathlib.Path(inputf)
    file_dir = input_path.parent
    allc_name = 'allc_' + input_path.name.split('.')[0] + 'tsv.gz'
    output_path = str(file_dir / allc_name)
    output_filehandler = gzip.open(output_path, 'wt')

    # process mpileup result
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    mc_sites = {'C', 'G'}
    context_len = num_upstr_bases + 1 + num_downstr_bases
    cur_chrom = ""
    line_counts = 0
    out = ""
    seq = None
    chr_out_pos_list = []
    cur_out_pos = 0
    for line in result_handle:
        fields = line.split("\t")
        # if chrom changed, read whole chrom seq from fasta
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            chr_out_pos_list.append((cur_chrom, str(cur_out_pos)))
            # get seq for cur_chrom
            seq = get_chromosome_sequence(reference_fasta, fai_df, cur_chrom)
            if seq is not None:
                seq = seq.upper()

        if seq is None:
            continue
        if fields[2] not in mc_sites:
            continue

        # indels
        read_bases = fields[4]
        incons_basecalls = read_bases.count("+") + read_bases.count("-")
        if incons_basecalls > 0:
            read_bases_no_indel = ""
            index = 0
            prev_index = 0
            while index < len(read_bases):
                if read_bases[index] == "+" or read_bases[index] == "-":
                    # get insert size
                    indel_size = ""
                    ind = index + 1
                    while True:
                        try:
                            int(read_bases[ind])
                            indel_size += read_bases[ind]
                            ind += 1
                        except:
                            break
                    try:
                        # sometimes +/- does not follow by a number and
                        # it should be ignored
                        indel_size = int(indel_size)
                    except:
                        index += 1
                        continue
                    read_bases_no_indel += read_bases[prev_index:index]
                    index = ind + indel_size
                    prev_index = index
                else:
                    index += 1
            read_bases_no_indel += read_bases[prev_index:index]
            fields[4] = read_bases_no_indel

        # count converted and unconverted bases
        if fields[2] == "C":
            pos = int(fields[1]) - 1
            try:
                context = seq[(pos - num_upstr_bases):(pos + num_downstr_bases + 1)]
            except:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(".")
            converted_c = fields[4].count("T")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = "\t".join([cur_chrom, str(pos + 1), "+", context,
                                  str(unconverted_c), str(cov), "1"]) + "\n"
                out += data
                cur_out_pos += len(data)

        elif fields[2] == "G":
            pos = int(fields[1]) - 1
            try:
                context = "".join([complement[base]
                                   for base in reversed(
                        seq[(pos - num_downstr_bases):(pos + num_upstr_bases + 1)]
                    )]
                                  )
            except:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(",")
            converted_c = fields[4].count("a")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = "\t".join([cur_chrom, str(pos + 1), "-", context,
                                  str(unconverted_c), str(cov), "1"]) + "\n"
                out += data
                cur_out_pos += len(data)

        if line_counts > buffer_line_number:
            output_filehandler.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_filehandler.write(out)
        line_counts = 0
        out = ""
    result_handle.close()
    output_filehandler.close()

    with open(output_path+'.idx', 'w') as idx_f:
        for (chrom, out_pos) in chr_out_pos_list:
            idx_f.write(f'{chrom}_{out_pos}\n')

    return 0
