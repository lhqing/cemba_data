from .bigwig import frag_to_bw_batch
from .dump_fragment_bed import dump_frags
from .macs2 import macs2


def atac_bulk_pipeline(cell_group_path,
                       output_dir_path,
                       sample_snap_path,
                       chrom_size_path,
                       species,
                       cpu=10,
                       remove_temp=True,
                       **macs2_kws):
    """
    Given a clustering assignment table and related SNAP files, do:
    1. extract fragment bed files per cluster
    2. generate BigWig file per cluster
    3. call peaks using MACS2 per cluster

    Parameters
    ----------
    cell_group_path:
        first row is header, names of cluster columns will be used as output names;
        Column 0 is sample name;
        Column 1 is cell barcode;
        Column 2, ..., n are level of cell group/clusters;
    output_dir_path
        Output directory, each cluster col will be a sub-dir
    sample_snap_path
        no header;
        Column 0 is sample name;
        Column 1 is sample SNAP file path;
    cpu
        Number of cpu to parallel
    chrom_size_path
        chromosome size file path
    species
        hs or mm
    remove_temp

    Returns
    -------

    """
    # extract cluster fragments
    frag_bed_path_list = dump_frags(
        cell_group_path,
        output_dir_path,
        sample_snap_path,
        cpu=cpu)

    # fragments to bigwig
    frag_to_bw_batch(
        frag_bed_path_list=frag_bed_path_list,
        chrom_size_path=chrom_size_path,
        remove_temp=remove_temp,
        cpu=cpu)

    # call peaks
    macs2(
        frag_bed_path_list=frag_bed_path_list,
        cpu=cpu,
        species=species,
        remove_temp=remove_temp,
        **macs2_kws)
