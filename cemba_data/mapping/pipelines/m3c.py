def m3c_config_str(config):
    """Change the dtype of parameters and make a appropriate string"""
    int_parameters = {
        'overlap': 6,
        'r1_left_cut': 10,
        'r1_right_cut': 10,
        'r2_left_cut': 10,
        'r2_right_cut': 10,
        'quality_threshold': 20,
        'length_threshold': 30,
        'total_read_pairs_min': 1,
        'total_read_pairs_max': 6000000,
        'mapq_threshold': 10,
        'num_upstr_bases': 0,
        'num_downstr_bases': 2,
        'compress_level': 5,
        'split_left_size': 40,
        'split_right_size': 40,
        'split_middle_min_size': 30,
        'min_gap': 1000,
        'trim_on_both_end': 5
    }

    str_parameters = {
        'mode': 'mc',
        'barcode_version': 'required',
        'r1_adapter': 'AGATCGGAAGAGCACACGTCTGAAC',
        'r2_adapter': 'AGATCGGAAGAGCGTCGTGTAGGGA',
        'bismark_reference': 'required',
        'reference_fasta': 'required',
        'mc_stat_feature': 'CHN CGN CCC',
        'mc_stat_alias': 'mCH mCG mCCC',
        'chrom_size_path': 'required'
    }

    typed_config = {}
    for k, default in int_parameters.items():
        if k in config:
            typed_config[k] = int(config[k])
        else:
            if default != 'required':
                typed_config[k] = default
            else:
                raise ValueError(f'Required parameter {k} not found in config. '
                                 f'You can print the newest mapping config template via "yap default-mapping-config".')

    for k, default in str_parameters.items():
        if k in config:
            typed_config[k] = f"'{config[k]}'"
        else:
            if default != 'required':
                typed_config[k] = f"'{default}'"
            else:
                raise ValueError(f'Required parameter {k} not found in config. '
                                 f'You can print the newest mapping config template via "yap default-mapping-config".')

    config_str = ""
    for k, v in typed_config.items():
        config_str += f"{k} = {v}\n"
    return config_str
