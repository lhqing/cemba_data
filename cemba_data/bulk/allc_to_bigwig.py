import pathlib

import pandas as pd


def generate_bigwig(output_dir_path,
                    chrom_size_path,
                    cpu,
                    mc_contexts_list,
                    bin_size):
    output_dir = pathlib.Path(output_dir_path).absolute()
    group_table_path = output_dir / 'GROUP_TABLE.csv'

    group_table = pd.read_csv(group_table_path, index_col=0)
    group_levels = group_table.columns[1:]
    qsub_dir = output_dir / 'qsub'

    bigwig_dir = output_dir / 'BW'
    bigwig_dir.mkdir(exist_ok=True)

    command_records = []
    for group_level in group_levels:
        _group_output_dir = bigwig_dir / group_level
        _group_output_dir.mkdir(exist_ok=True)

        file_records = []
        group_dir = output_dir / group_level
        group_allc_records = group_dir / 'ALLC_table.csv'
        group_allc_df = pd.read_csv(group_allc_records, index_col=0)

        for allc_path in group_allc_df['AllcPath']:
            mc_contexts_str = ' '.join(mc_contexts_list)
            output_prefix = _group_output_dir / pathlib.Path(allc_path).name[:-12]
            output_path = ''  # TODO

            command = f'allcools allc-to-bigwig ' \
                      f'--allc_path {allc_path} ' \
                      f'--output_prefix {output_prefix} ' \
                      f'--chrom_size_path {chrom_size_path} ' \
                      f'--mc_contexts {mc_contexts_str} ' \
                      f'--bin_size {bin_size} ' \
                      f'--remove_additional_chrom ' \
                      f'--cpu {cpu}'
            command_records.append(command)
            file_records.append(output_path)

        # TODO group_allc_records will be updated
        # group_allc_df[''] = file_records
        # group_allc_df.to_csv(group_allc_records)

    with open(qsub_dir / f'_bigwig.commands.txt', 'w') as f:
        f.write('\n'.join(command_records))
    return
