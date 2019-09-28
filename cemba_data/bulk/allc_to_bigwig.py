import json

def generate_bigwig(output_dir,
                    sub_dir_name,
                    chrom_size_path,
                    cpu,
                    mc_contexts_list,
                    bin_size):
    qsub_dir = output_dir / 'qsub/bigwig'
    qsub_dir.mkdir(exist_ok=True)

    _output_dir = output_dir / sub_dir_name
    output_records_path = _output_dir / f'merge_{sub_dir_name}.records.json'
    with open(output_records_path) as f:
        allc_dict = json.load(f)
    command_records = []
    output_records = {}
    for cluster, allc_path in allc_dict.items():
        output_prefix = _output_dir / cluster

        if allc_path.is_symlink():
            real_allc = allc_path.resolve()
            bw_files = real_allc.parent.glob(f'{cluster}*bw')
            for bw_file_path in bw_files:
                bw_suffix = '.'.join(bw_file_path.name.split('.')[-3:])
                command = f'ln -s {bw_file_path} {output_prefix}.{bw_suffix}'
                command_records.append(command)
        else:
            mc_contexts_str = ' '.join(mc_contexts_list)
            command = f'allcools allc-to-bigwig ' \
                      f'--allc_path {allc_path} ' \
                      f'--output_prefix {output_prefix} ' \
                      f'--chrom_size_path {chrom_size_path} ' \
                      f'--mc_contexts {mc_contexts_str} ' \
                      f'--bin_size {bin_size} ' \
                      f'--remove_additional_chrom ' \
                      f'--cpu {cpu}'
            command_records.append(command)
    with open(qsub_dir / f'bigwig_{sub_dir_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(_output_dir / f'bigwig_{sub_dir_name}.records.json', 'w') as f:
        json.dump(output_records, f)
