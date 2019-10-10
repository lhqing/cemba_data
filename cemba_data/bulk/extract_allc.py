import pathlib

import pandas as pd


def extract_strand_merged_cg(output_dir_path, chrom_size_path, mc_context='CGN', cpu=5):
    output_dir = pathlib.Path(output_dir_path).absolute()
    group_table_path = output_dir / 'GROUP_TABLE.csv'

    group_table = pd.read_csv(group_table_path, index_col=0)
    group_levels = group_table.columns[1:]
    qsub_dir = output_dir / 'qsub'

    commands = []
    for group_level in group_levels:
        file_records = []
        group_dir = output_dir / group_level
        group_allc_records = group_dir / 'ALLC_table.csv'
        group_allc_df = pd.read_csv(group_allc_records, index_col=0)

        for allc_path in group_allc_df['AllcPath']:
            output_prefix = str(allc_path[:-12])
            output_path = output_prefix + '.CGN-Merge.allc.tsv.gz'
            cmd = f'allcools extract-allc ' \
                  f'--allc_path {allc_path} ' \
                  f'--output_prefix {output_prefix} ' \
                  f'--mc_contexts {mc_context} ' \
                  f'--chrom_size_path {chrom_size_path} ' \
                  f'--strandness merge ' \
                  f'--output_format allc ' \
                  f'--cov_cutoff 99999 ' \
                  f'--cpu {cpu}'
            commands.append(cmd)
            file_records.append(output_path)
        # group_allc_records will be updated
        group_allc_df['mCGExtract'] = file_records
        group_allc_df.to_csv(group_allc_records)

    command_file_path = qsub_dir / '_extract_cg_commands.txt'
    with open(command_file_path, 'w') as f:
        f.write('\n'.join(commands))
    return command_file_path
