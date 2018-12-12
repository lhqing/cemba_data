import pathlib
import os
import json
import glob
import pandas as pd
import argparse
from ..mapping.pipeline import validate_fastq_dataframe
import logging

# logger
log = logging.getLogger()


def get_fastq_dataframe(file_path, name_pattern,
                        uid_field='uid', lane_field='lane', read_type_field='read-type',
                        out_dir=None):
    """
    Generate fastq_dataframe for pipeline input.

    Parameters
    ----------
    file_path
        Accept 1. path pattern, 2. path list, 3. path of one file contain all the paths.
    name_pattern
        Name pattern of file paths, format example: "Field1_Field2_Field3_..._Fieldn".
        Both file name (after remove suffix part after frist '.') and
        name pattern are split by '_' and match together.
        If the length after split is not the same, error will be raised.
    uid_field
        field(s) in name pattern that can combine to uid for multiplexed fileset. will be linked by '-'
    lane_field
        field in name pattern that stand for lane
    read_type_field
        field in name pattern that stand for read type (R1, R2)
    out_dir
        output dir path for fastq dataframe

    Returns
    -------
        fastq_dataframe for pipeline input.
    """
    if isinstance(file_path, str) and ('*' in file_path):
        file_path = [str(pathlib.Path(p).absolute()) for p in glob.glob(file_path)]
    elif isinstance(file_path, list):
        pass
    else:
        with open(file_path) as f:
            file_path = [l.strip() for l in f]

    name_patterns = name_pattern.lower().split('_')

    uid_fields = uid_field.lower().split('_')
    lane_fields = lane_field.lower().split('_')
    read_type_fields = read_type_field.lower().split('_')
    for field in (uid_fields + lane_fields + read_type_fields):
        if field not in name_pattern:
            raise ValueError(f'field name {field} not in name pattern {name_pattern}')

    fastq_data = []
    broken_names = []
    for fastq in file_path:
        fastq_name_list = fastq.split('/')[-1].split('.')[0].split('_')
        if len(fastq_name_list) != len(name_patterns):
            broken_names.append(fastq)
            continue
        name_dict = {p: v for p, v in zip(name_patterns, fastq_name_list)}
        name_dict['fastq_path'] = fastq
        name_series = pd.Series(name_dict)
        fastq_data.append(name_series)

    fastq_df = pd.DataFrame(fastq_data)
    log.info(len(broken_names), 'broken names.')
    log.info(fastq_df.shape[0], 'valid fastq names.')
    if fastq_df.shape[0] == 0:
        log.info('No fastq name remained, check if the name pattern is correct.')
        return None

    fastq_df['uid'] = fastq_df[uid_fields].apply(lambda i: '-'.join(i.tolist()), axis=1)
    fastq_df['lane'] = fastq_df[lane_fields].apply(lambda i: '-'.join(i.tolist()), axis=1)
    fastq_df['read_type'] = fastq_df[read_type_fields].apply(lambda i: '-'.join(i.tolist()), axis=1)
    # make sure UID is unique
    for _, df in fastq_df.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError(f'UID based on definition "{uid_field}" are not unique.')

    for remove_col in ['*', 'read-type']:
        if remove_col in fastq_df.columns:
            del fastq_df[remove_col]

    if out_dir is not None:
        fastq_df.to_csv(pathlib.Path(out_dir) / 'fastq_dataframe.tsv.gz',
                        index=None, sep='\t', compression='gzip')
        if len(broken_names) != 0:
            with open(pathlib.Path(out_dir) / 'fastq_dataframe.broken_names.txt', 'w') as f:
                f.write('\n'.join(broken_names))
    return fastq_df


def batch_pipeline(fastq_dataframe, out_dir, config_path):
    """
    Parallel version of mapping pipeline.
    Only generate a command.json file, not actually running the mapping.
    Use yap qsub to submit the mapping jobs.

    Parameters
    ----------
    fastq_dataframe
        Input dataframe for all fastq file path and metadata.
        Must include columns: uid, read_type, index_name, lane
    out_dir
        pipeline universal out_dir
    config_path
        pipeline universal config
    """

    # prepare out_dir
    out_dir = pathlib.Path(out_dir).absolute()
    out_dir.mkdir(parents=True)
    # prepare config
    if config_path is None:
        config_path = os.path.dirname(__file__) + '/mapping_config.ini'
    with open(out_dir / 'mapping_config.ini', 'w') as wf, open(config_path, 'r') as f:
        for line in f:
            wf.write(line)
    config_path = out_dir / 'mapping_config.ini'
    # prepare fastq dataframe
    fastq_dataframe = validate_fastq_dataframe(fastq_dataframe)
    fastq_dataframe_path = out_dir / 'fastq_dataframe.tsv.gz'
    fastq_dataframe.to_csv(fastq_dataframe_path,
                           sep='\t', compression='gzip', index=None)

    # generate command json
    cmd_list = []
    for uid, sub_df in fastq_dataframe.groupby('uid'):
        uid_out_dir = out_dir / uid
        uid_out_dir.mkdir()
        uid_fastq_dataframe_path = uid_out_dir / 'fastq_dataframe.tsv.gz'
        sub_df.to_csv(uid_fastq_dataframe_path,
                      sep='\t', compression='gzip', index=None)
        uid_cmd = f'yap mapping --fastq_dataframe {uid_fastq_dataframe_path} ' \
                  f'--out_dir {uid_out_dir} --config_path {config_path}'
        command_dict = {
            'command': uid_cmd,
            '-pe smp': 20,  # cpu for each command
            # TODO set cup number more clever
        }
        cmd_list.append(command_dict)
    cmd_json_path = out_dir / 'command.json'
    with open(cmd_json_path, 'w') as f:
        json.dump(cmd_list, f)

    qsub_command = f'yap qsub --working_dir {out_dir} ' \
                   f'--project_name mapping ' \
                   f'--command_file_path {cmd_json_path} ' \
                   f'--total_cpu 160 --submission_gap 2 --qstat_gap 60'

    log.info(f"""
    The output directory is prepared. All related files have been copied to that directory.
    ---------------------------------------------------------------------------------------
    - Output directory: {out_dir}
    - Configuration file: {config_path}
    - FASTQ dataframe: {fastq_dataframe_path}
    - Yap mapping command list: {cmd_json_path}
    
    To run the command list using qsub, use "yap qsub" like this:
    
    {qsub_command}
    
    Modify the qsub parameters if you need. See "yap qsub -h" for help.
    """)

    return str(cmd_json_path)


def batch_pipeline_register_subparser(subparser):
    parser = subparser.add_parser('mapping-qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="Mapping pipeline from multiplexed FASTQ file to ALLC file.")
    parser.set_defaults(func=batch_pipeline)

    parser_req = parser.add_argument_group("Required inputs")

    parser_req.add_argument(
        "--fastq_dataframe",
        type=str,
        required=True,
        help="Path of fastq dataframe, can be generate with yap fastq_dataframe"
    )

    parser_req.add_argument(
        "--out_dir",
        type=str,
        required=True,
        help="Pipeline output directory, if not exist, will create recursively."
    )

    parser_req.add_argument(
        "--config_path",
        type=str,
        required=True,
        default=None,
        help="Pipeline configuration (.ini) file path. "
             "You can use 'yap default-mapping-config' to print out default config can modify it."
    )
    return
