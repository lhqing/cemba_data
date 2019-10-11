import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pybedtools


def fragment_to_bigwig(fragment_bed_path,
                       chrom_size_path,
                       output_prefix,
                       scale_factor=1e6,
                       remove_bg=True,
                       sort_mem='2%',
                       sort_cpu=1):
    """Cluster fragment bed file to bigwig file"""
    output_prefix = str(output_prefix).rstrip('.')
    out_bw_path = output_prefix + '.bw'

    print('Sort fragment')
    out_sort_path = str(out_bw_path) + 'temp.sort.bed'
    try:
        subprocess.run(f'zcat {fragment_bed_path} | sort -k1,1 -k2,2n '
                       f'-S {sort_mem} --parallel {int(sort_cpu)} > {out_sort_path}',
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8',
                       check=True, shell=True)
        # count reads
        p = subprocess.run(['wc', '-l', out_sort_path],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8',
                           check=True)
        total_reads = int(p.stdout.split(' ')[0])
        scale = scale_factor / total_reads
    except subprocess.CalledProcessError as e:
        print(e.stderr)
        raise e

    print('Calculate COV')
    frag_sorted_bed = pybedtools.BedTool(out_sort_path)
    cov_bed = frag_sorted_bed.genome_coverage(scale=scale, bg=True, g=chrom_size_path)
    out_bg_path = output_prefix + '.bg'
    cov_bed.saveas(str(out_bg_path))

    print('Generate bigwig')
    subprocess.run(['bedGraphToBigWig', out_bg_path, chrom_size_path, out_bw_path],
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8', check=True)
    if remove_bg:
        subprocess.run(['rm', '-f', out_bg_path, out_sort_path], check=True)
    else:
        subprocess.run(['mv', '-f', out_sort_path], check=True)
    return out_bw_path


def frag_to_bw_batch(frag_bed_path_list, chrom_size_path, remove_temp=True, cpu=1):
    future_list = []
    with ProcessPoolExecutor(cpu) as executor:
        for fragment_bed_path in frag_bed_path_list:
            future = executor.submit(fragment_to_bigwig,
                                     fragment_bed_path=fragment_bed_path,
                                     chrom_size_path=chrom_size_path,
                                     output_prefix=str(fragment_bed_path)[:-7],
                                     scale_factor=1e6,
                                     remove_bg=remove_temp,
                                     sort_mem='2%',
                                     sort_cpu=1)
            future_list.append(future)
    for future in as_completed(future_list):
        future.result()
    return
