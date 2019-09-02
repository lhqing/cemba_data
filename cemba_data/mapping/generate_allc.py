import pathlib
from .utilities import get_configuration
import pandas as pd


def generate_allc(output_dir, config):
    output_dir = pathlib.Path(output_dir)
    config = get_configuration(config)

    bam_records = pd.read_csv(output_dir / 'final_bam.records.csv',
                              index_col=['uid', 'index_name', 'read_type'],
                              squeeze=True)
    reference_fasta = config['callMethylation']['reference_fasta']
    num_upstr_bases = config['callMethylation']['num_upstr_bases']
    num_downstr_bases = config['callMethylation']['num_downstr_bases']
    compress_level = config['callMethylation']['compress_level']
    min_mapq = config['callMethylation']['min_mapq']
    min_base_quality = config['callMethylation']['min_base_quality']

    records = []
    command_list = []
    for (uid, index_name, read_type), bismark_bam_path in bam_records.iteritems():
        # file path
        output_allc_path = output_dir / f'{uid}-{index_name}.allc.tsv.gz'
        # command
        command = f'allcools bam-to-allc ' \
                  f'--bam_path {bismark_bam_path} ' \
                  f'--reference_fasta {reference_fasta} ' \
                  f'--output_path {output_allc_path} ' \
                  f'--cpu 1 ' \
                  f'--num_upstr_bases {num_upstr_bases} ' \
                  f'--num_downstr_bases {num_downstr_bases} ' \
                  f'--min_mapq {min_mapq} ' \
                  f'--min_base_quality {min_base_quality} ' \
                  f'--compress_level {compress_level} ' \
                  f'--save_count_df'
        records.append([uid, index_name, read_type, output_allc_path])
        command_list.append(command)

    with open(output_dir / 'generate_allc.command.txt', 'w') as f:
        f.write('\n'.join(command_list))
    record_df = pd.DataFrame(records,
                             columns=['uid', 'index_name', 'read_type', 'allc_path'])
    record_df.to_csv(output_dir / 'generate_allc.records.csv', index=None)
    return record_df, command_list
