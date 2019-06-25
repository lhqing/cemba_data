import pathlib
import shlex
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def runner(cmd):
    print(cmd)
    p = subprocess.run(shlex.split(cmd), check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       encoding='utf8')
    return p


def merge_cluster_bam(bam_series, cluster_series, output_dir, cpu=1):
    """
    Merge single cell bam to cluster level

    Parameters
    ----------
    bam_series
        index is cell_id, value is bam path
    cluster_series
        index is cell_id, value is cluster assignment
    output_dir
        Output Directory
    cpu
        number of CPU to use

    Returns
    -------

    """
    output_dir = pathlib.Path(output_dir).absolute()
    with ProcessPoolExecutor(cpu) as executor:
        futures = []
        for cluster, sub_list in bam_series.groupby(cluster_series):
            cell_list_path = output_dir / f'{cluster}.cell_bam_list'
            with open(cell_list_path, 'w') as f:
                for path in sub_list:
                    f.write(path + '\n')
            out_path = output_dir / f'{cluster}.merge.bam'
            cmd = f'samtools merge -b {cell_list_path} -@ 2 {out_path}'
            future = executor.submit(runner, cmd)
            futures.append(future)

        for future in as_completed(futures):
            try:
                print(future.result())
            except Exception as e:
                raise e
    return
