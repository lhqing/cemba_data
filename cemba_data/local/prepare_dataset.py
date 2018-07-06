import configparser
import os

ref_path_config = configparser.ConfigParser()
ref_path_config.read(os.path.dirname(__file__) + '/config_ref_path.ini')


def batch_map_to_region(allc_files, out_dir,
                        region_bed_path, region_name, genome_size_path,
                        context_pattern, max_cov_cutoff,
                        remove_tmp=True, tmp_compression=False):
    return


def assembl_dataset():
    return


def prepare_dataset():
    # cell df and filter parms (defalut none if df is whole dataset)
    # 1. use batch_map_to_region get all files
    # 2. assemble_dataset to h5 file

    return

