"""
A sbatch wrapper for stampede2
See stampede2 doc:
https://portal.tacc.utexas.edu/user-guides/stampede2#running-sbatch
"""

import json
import pathlib
import re
import shlex
import subprocess
import time

import pandas as pd

import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])

# see stampede2 doc https://portal.tacc.utexas.edu/user-guides/stampede2#running-queues
# name: max_jobs
STAMPEDE2_QUEUES = {
    'development': 1,
    'normal': 50,
    'large': 5,
    'long': 2,
    'flat-quadrant': 5,
    'skx-dev': 1,
    'skx-normal': 25,
    'skx-large': 3
}


def get_job_id(sbatch_result):
    """
    parse the sbatch output, check status, and get job id

    Parameters
    ----------
    sbatch_result

    Returns
    -------
    sbatch job_id
    """
    job_id = None
    for line in sbatch_result:
        line = line.strip()
        if line.startswith('-->'):
            # status line
            if line.endswith("OK"):
                continue
            else:
                print(sbatch_result)
                raise ValueError("sbatch output is abnormal, see information above")
        elif line.startswith('Submitted batch job '):
            # job id line
            job_id = line.split(' ')[-1]
        else:
            pass
    if job_id is None:
        print(sbatch_result)
        raise ValueError('Can not get job id from sbatch output, see information above')
    return job_id


def submit_sbatch(job_script_path):
    """
    submit sbatch job and return job id

    Parameters
    ----------
    job_script_path

    Returns
    -------

    """
    p = subprocess.run(['sbatch', job_script_path],
                       check=True,
                       encoding='utf8',
                       stdout=subprocess.PIPE)
    job_id = get_job_id(p.stdout)
    print(f'Submit job script: {job_script_path}. Job ID is {job_id}.')
    return job_id


def squeue():
    """
    check current running job

    Returns
    -------
    squeue results in a pd.DataFrame
    """
    user_name = subprocess.run(['whoami'],
                               check=True,
                               encoding='utf8',
                               stdout=subprocess.PIPE).stdout.strip()
    squeue_result = subprocess.run(['squeue', '-u', user_name],
                                   check=True,
                                   encoding='utf8',
                                   stdout=subprocess.PIPE).stdout

    sep_pattern = re.compile(r' +')
    contents = [sep_pattern.split(line.strip()) for line in squeue_result]
    squeue_df = pd.DataFrame(contents[1:],
                             columns=contents[0]).set_index('JOBID')
    print('Current squeue output:')
    print(squeue_result, end='\n\n')
    return squeue_df


def check_targets(targets):
    """
    Check whether all target file paths exist, return bool

    Parameters
    ----------
    targets

    Returns
    -------
    True or False
    """
    for path in targets:
        if not pathlib.Path(path).exists():
            return False
    return True


def make_sbatch_script_files(commands, sbatch_dir, name_prefix, queue, time_str, email, email_type):
    """See stampede2 doc: https://portal.tacc.utexas.edu/user-guides/stampede2#running-sbatch"""
    with open(PACKAGE_DIR / 'files/sbatch_template.txt') as f:
        sbatch_template = f.read()
    sbatch_dir = pathlib.Path(sbatch_dir)

    if email is not None:
        email_str = f'#SBATCH --mail-user={email}'
        email_type_str = f'#SBATCH --mail-type={email_type}'
    else:
        email_str = ''
        email_type_str = ''

    command_to_job_script_paths = {}
    for i, command in enumerate(commands):
        job_name = f'{name_prefix}_{i}'
        sbatch_script = sbatch_template.format(
            job_name=job_name,
            queue=queue,
            time_str=time_str,
            email_str=email_str,
            email_type_str=email_type_str,
            command=command
        )
        job_script_path = sbatch_dir / f'{job_name}.sh'
        with open(job_script_path, 'w') as f:
            f.write(sbatch_script)
        command_to_job_script_paths[command] = job_script_path
    return command_to_job_script_paths


def sacct(jobs):
    sep_pattern = re.compile(r' +')
    sacct_cmd = f'sacct -j {",".join(jobs)} ' \
                f'--format=jobid,jobname,partition,alloccpus,elapsed,state,exitcode'
    sacct_result = subprocess.run(shlex.split(sacct_cmd),
                                  check=True,
                                  encoding='utf8',
                                  stdout=subprocess.PIPE).stdout
    lines = []
    header = ''
    col_starts = []
    for i, line in enumerate(sacct_result):
        line = line.rstrip('\n')
        if i == 0:
            header = line
        elif i == 1:
            # the second row indicates col width, use it to determine col_starts
            col_width = [len(s) for s in sep_pattern.split(line)]
            cur_pos = 0
            col_starts = [0]
            for length in col_width:
                cur_pos += (length + 1)
                col_starts.append(cur_pos)
        else:
            lines.append(line)

    columns = [header[col_starts[i]:col_starts[i + 1]].strip() for i in range(len(col_starts) - 1)]
    data = []
    for line in lines:
        ll = [line[col_starts[i]:col_starts[i + 1]].strip() for i in range(len(col_starts) - 1)]
        data.append(ll)
    sacct_data = pd.DataFrame(data, columns=columns).set_index('JobID')
    sacct_data = sacct_data[~sacct_data.index.str.endswith('bat+')].copy()
    sacct_data['Success'] = sacct_data['ExitCode'] == '0:0'
    return sacct_data


def sbatch_submitter(name, commands_to_targets, output_dir, time_str, queue='skx-normal',
                     email=None, email_type='fail', max_jobs=None):
    with open(commands_to_targets) as f:
        commands_to_targets = json.load(commands_to_targets)
    # set name
    name = name.replace(' ', '_')

    # check queue
    queue = queue.lower()
    if queue not in STAMPEDE2_QUEUES:
        raise KeyError(f'queue name {queue} not found in STAMPEDE2_QUEUES, '
                       f'available queues are {list(STAMPEDE2_QUEUES.keys())}')

    # set max_jobs
    _max_jobs = STAMPEDE2_QUEUES[queue]
    if max_jobs is None:
        max_jobs = _max_jobs
    else:
        max_jobs = min(max_jobs, _max_jobs)
    print(f'Max concurrent sbatch jobs {max_jobs}, stampede2 allows {STAMPEDE2_QUEUES[queue]}.')

    # create job script files
    sbatch_dir = pathlib.Path(output_dir) / f'{name}_sbatch'
    sbatch_dir.mkdir(exist_ok=True)
    command_to_job_script_paths = make_sbatch_script_files(
        commands=commands_to_targets.keys(),
        sbatch_dir=sbatch_dir,
        name_prefix=name,
        queue=queue,
        time_str=time_str,
        email=email,
        email_type=email_type)

    # prepare submission
    running_job_id_set = set()  # sbatch_id
    finished_job_and_targets = {}  # {job_id: target status}
    job_script_path_to_targets = {script_path: commands_to_targets[command]  # {job_script_path: target_paths}
                                  for command, script_path in command_to_job_script_paths.items()}
    queue_job_path_list = list(job_script_path_to_targets.keys()).copy()
    job_id_to_script_path = {}

    # start submission
    sleepy = 30
    flag_path = sbatch_dir / 'RUNNING_SIGNAL'
    with open(flag_path, 'w') as f:
        f.write('')
    while (len(queue_job_path_list) != 0) or (len(running_job_id_set) != 0):
        if not flag_path.exists():
            # break if flag missing
            break

        # squeue and update running job status
        squeue_df = squeue()
        remaining_slots = max_jobs - squeue_df.shape[0]
        # the max_jobs is apply to user level, not to the current submitter level
        if remaining_slots > 0:
            # things are getting done, weak up
            sleepy = 30

            # check running jobs
            new_running_job_id_set = set()
            for job_id in running_job_id_set:
                if job_id not in squeue_df.index:
                    # running job finished
                    script_path = job_id_to_script_path[job_id]
                    finished_job_and_targets[job_id] = check_targets(job_script_path_to_targets[script_path])
                    # status will be judged in the end
                else:
                    # still running
                    new_running_job_id_set.add(job_id)
                    pass
            running_job_id_set = new_running_job_id_set
            print(f'{len(running_job_id_set)} running job IDs: {", ".join(running_job_id_set)}')

            # submit new jobs
            print(f'Going to submit {min(len(queue_job_path_list), remaining_slots)} new jobs')
            while (remaining_slots > 0) and (len(queue_job_path_list) > 0):
                script_path = queue_job_path_list.pop()
                # skip if target already exist
                targets = job_script_path_to_targets[script_path]
                if check_targets(targets):
                    print(f'Target file exist for script {script_path}')
                    continue
                else:
                    job_id = submit_sbatch(script_path)
                    running_job_id_set.add(job_id)
                    job_id_to_script_path[job_id] = script_path
                    remaining_slots -= 1

        # sleep
        sleepy += 30
        sleepy = min(300, sleepy)
        time.sleep(sleepy)

    # only check status if something has finished
    if len(finished_job_and_targets) > 0:
        # check status
        chunk_size = 100
        finished_job_ids = list(finished_job_and_targets.keys())
        stats = []
        for i in range(0, len(finished_job_ids), chunk_size):
            job_chunk = finished_job_ids[i: i + chunk_size]
            stats.append(sacct(job_chunk))
        sacct_df = pd.concat(stats)
        sacct_df['TargetExist'] = sacct_df.index.map(finished_job_and_targets)
        sacct_df['ScriptPath'] = sacct_df.index.map(job_id_to_script_path)
        sacct_df['Success'] = sacct_df['Success'] & sacct_df['TargetExist']  # recheck success

        # append old status if exist
        sacct_path = sbatch_dir / 'sacct.csv.gz'
        if sacct_path.exists():
            df = pd.read_csv(sacct_path, index_col=0)
            # remove resubmitted jobs, so each ScriptPath only has one row
            df = df[~df['ScriptPath'].isin(set(sacct_df['ScriptPath']))]
            sacct_df = pd.concat([df, sacct_df])
        sacct_df.to_csv(sacct_path)

        # print time elapsed
        times = pd.to_timedelta(sacct_df[sacct_df['Success']]['Elapsed'])
        print('For all the successful jobs')
        print('The minimum job elapsed:', times.min())
        print('The average job elapsed:', times.mean())
        print('The maximum job elapsed:', times.max())

        # print success
        print(f"{sacct_df['Success'].sum()} / {sacct_df.shape[0]} succeeded")

    # delete flag
    subprocess.run(shlex.split(f'rm -f {flag_path}'), check=True)
    return
