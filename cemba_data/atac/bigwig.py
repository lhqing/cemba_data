import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import pandas as pd
import pybedtools


def fragment_to_bigwig(fragment_bed_path,
                       chrom_size_path,
                       output_prefix,
                       scale_factor=1e6,
                       remove_bg=True,
                       sort_mem='2%',
                       sort_cpu=1):
    """Cluster fragment bed file to bigwig file"""
    # TODO this function has memory problem

    output_prefix = str(output_prefix).rstrip('.')
    out_bw_path = output_prefix + '.bw'
    if pathlib.Path(out_bw_path).exists():
        return

    frag_bed_df = pd.read_csv(fragment_bed_path, sep='\t', header=None,
                              names=['chrom', 'start', 'end', 'fragment'])
    frag_bed = pybedtools.BedTool.from_dataframe(frag_bed_df)
    frag_sorted_bed = frag_bed.sort(g=chrom_size_path)
    scale = scale_factor / frag_bed_df.shape[0]

    print('Calculate COV')
    cov_bed = frag_sorted_bed.genome_coverage(scale=scale, bg=True, g=chrom_size_path)
    out_path = output_prefix + '.bg'
    cov_bed.saveas(str(out_path))

    print('Generate bigwig')
    out_sort_path = str(out_path) + '.sort'
    with open(str(out_path) + '.sort', 'wb') as f:
        p = subprocess.run(['sort', '-k1,1', '-k2,2n',
                            '-S', sort_mem, '--parallel', str(sort_cpu), out_path],
                           stdout=subprocess.PIPE)
        f.write(p.stdout)
    subprocess.run(['bedGraphToBigWig', out_sort_path, chrom_size_path, out_bw_path], check=True)
    if remove_bg:
        subprocess.run(['rm', '-f', out_path, out_sort_path], check=True)
    else:
        subprocess.run(['mv', '-f', out_sort_path, out_path], check=True)
    return out_bw_path


def frag_to_bw_batch(cell_group_path, output_dir_path, chrom_size_path, cpu=1):
    cell_group_table = pd.read_csv(cell_group_path, index_col=None)
    output_dir = pathlib.Path(output_dir_path).absolute()
    future_list = []
    with ProcessPoolExecutor(cpu) as executor:
        for col in cell_group_table.columns[2:]:
            col_dir = output_dir / col
            for fragment_bed_path in col_dir.glob('*.bed.gz'):
                future = executor.submit(fragment_to_bigwig,
                                         fragment_bed_path=fragment_bed_path,
                                         chrom_size_path=chrom_size_path,
                                         output_prefix=str(fragment_bed_path)[:-7],
                                         scale_factor=1e6,
                                         remove_bg=True,
                                         sort_mem='2%',
                                         sort_cpu=1)
                future_list.append(future)
    for future in as_completed(future_list):
        future.result()
    return
