from .dump_fragment_bed import dump_frags
from .bigwig import frag_to_bw_batch


def atac_bulk_pipeline(cell_group_path, output_dir_path, sample_snap_path, chrom_size_path, cpu):
    dump_frags(cell_group_path, output_dir_path, sample_snap_path, cpu=cpu)

    frag_to_bw_batch(cell_group_path, output_dir_path, chrom_size_path, cpu=cpu)

    # call peaks
