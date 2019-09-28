
def extract_cg(output_dir, sub_dir_name, chrom_size_path, cpu, mc_contexts='CGN'):
    qsub_dir = output_dir / 'qsub/extract_cg'
    qsub_dir.mkdir(exist_ok=True)

    _output_dir = output_dir / sub_dir_name
    output_records_path = _output_dir / f'merge_{sub_dir_name}.records.json'
    with open(output_records_path) as f:
        allc_dict = json.load(f)
    command_records = []
    output_records = {}
    for cluster, allc_path in allc_dict.items():
        output_prefix = _output_dir / cluster
        output_path = _output_dir / (cluster + f'.{mc_contexts}-Merge.allc.tsv.gz')
        output_records[cluster] = output_path
        if allc_path.is_symlink():
            # for these smaller files, use real command and actual computation is clearer
            allc_path = allc_path.resolve()
        command = f'allcools extract-allc ' \
                  f'--allc_path {allc_path} ' \
                  f'--output_prefix {output_prefix} ' \
                  f'--mc_contexts {mc_contexts} ' \
                  f'--chrom_size_path {chrom_size_path} ' \
                  f'--strandness merge ' \
                  f'--output_format allc ' \
                  f'--cpu {cpu}'
        command_records.append(command)
    with open(qsub_dir / f'extract_cg_{sub_dir_name}.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    with open(_output_dir / f'extract_cg_{sub_dir_name}.records.json', 'w') as f:
        json.dump(output_records, f)

