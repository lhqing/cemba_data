import pathlib
import os
import json
import glob
import pandas as pd
from ...mapping.utilities import validate_fastq_dataframe
import logging

# logger
log = logging.getLogger()


def get_fastq_dataframe(file_path, output_path=None, skip_broken_name=False):
    """
    Generate fastq_dataframe for pipeline input.

    Parameters
    ----------
    file_path
        Accept 1. path pattern contain wildcard, 2. path list, 3. path of one file contain all the paths.
    output_path
        output path of the fastq dataframe
    skip_broken_name
        If true, ignore any unrecognized file names in file_path
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
    log.info(f'{len(file_path)} fastq file paths in input')

    fastq_data = []
    broken_names = []
    for path in file_path:
        path = pathlib.Path(path)
        try:
            *_, plate1, plate2, multi_field = path.name.split('-')
            plate_pos, _, lane, read_type, _ = multi_field.split('_')
            try:
                assert plate_pos[0] in 'ABCDEFGH'
                assert int(plate_pos[1:]) in list(range(1, 13))
                assert lane in {'L001', 'L002', 'L003', 'L004'}
                assert read_type in {'R1', 'R2'}
                assert plate1 != plate2
            except AssertionError:
                raise ValueError
        except ValueError:
            broken_names.append(path)
            if skip_broken_name:
                continue
            raise ValueError(f'Found unknown name pattern in path {path}')
        name_dict = dict(plate1=plate1,
                         plate2=plate2,
                         plate_pos=plate_pos,
                         lane=lane,
                         read_type=read_type,
                         fastq_path=path,
                         uid=f'{plate1}-{plate2}-{plate_pos}')
        name_series = pd.Series(name_dict)
        fastq_data.append(name_series)

    fastq_df = pd.DataFrame(fastq_data)
    log.info(f'{len(broken_names)} broken names.')
    log.info(f'{fastq_df.shape[0]} valid fastq names.')
    if fastq_df.shape[0] == 0:
        log.info('No fastq name remained, check if the name pattern is correct.')
        return None

    # make sure UID is unique
    for _, df in fastq_df.groupby(['lane', 'read_type']):
        if df['uid'].unique().size != df['uid'].size:
            raise ValueError(f'UID column is not unique.')
    if output_path is not None:
        if not output_path.endswith('tsv.gz'):
            output_path += 'tsv.gz'
        fastq_df.to_csv(output_path, index=None, sep='\t', compression='gzip')
        return
    else:
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
        pipeline universal output_dir
    config_path
        pipeline universal config
    """

    # prepare output_dir
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
                  f'--output_dir {uid_out_dir} --config_path {config_path}'
        command_dict = {
            'command': uid_cmd,
            'pe smp': 20,  # cpu for each command
            'l h_vmem': '5G'
        }
        cmd_list.append(command_dict)
    cmd_json_path = out_dir / 'command.json'
    with open(cmd_json_path, 'w') as f:
        json.dump(cmd_list, f)

    qsub_command = f'yap qsub --working_dir {out_dir} ' \
                   f'--project_name allc ' \
                   f'--command_file_path {cmd_json_path} ' \
                   f'--total_cpu 160'

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
