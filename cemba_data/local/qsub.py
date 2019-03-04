"""
SGE qsub system auto-submitter
# TODO add runtime parameter change, e.g. submitted using 100 CPU, change to 200 at runtime
# TODO add safe_terminate, e.g. terminate after all current running job have been finished, but not submit new job
# TODO deal with failed and incomplete json in resubmit
# TODO add option to profile (qacct) several jobs as reference for next submission
# TODO add test_run parameter to run N command and profile each and give report
"""

from subprocess import run, PIPE, CalledProcessError
import os
import json
import datetime
import re
import sys
import time
import shutil
import pandas as pd

SUBMISSION_GAP = 2
QSTAT_GAP = 30


def _default_command_dict(name, error_path, output_path, working_dir,
                          cpu=1, memory='2G', h_rt='99:99:99', s_rt='99:99:99'):
    """Default qsub command dict"""
    command_dict = {'command': None,
                    'N': name,
                    'V': '',
                    'pe smp': cpu,
                    'l h_vmem': memory,
                    'l h_rt': h_rt,
                    'l s_rt': s_rt,
                    'wd': working_dir,
                    'e': error_path,
                    'o': output_path}
    return command_dict


def _get_running_job_id_qstat(user_name, id_set):
    """
    Check qstat of user_name, given all submitted id_set

    Parameters
    ----------
    user_name
    id_set

    Returns
    -------
    Id list of all current running id
    """

    while True:
        tried = 0
        try:
            qstat_result_string = run(['qstat', '-u', user_name],
                                      stdout=PIPE, encoding='utf8', check=True).stdout
            break
        except CalledProcessError as e:
            # TODO: there seems to be a bug for qstat occasionally happen in mapping,
            # qstat return through an unknown error about
            # ERROR: failed receiving gdi request response for mid=1 (got syncron message receive timeout error)
            tried += 1
            time.sleep(30)
            if tried > 10:
                print(f'Qstat failed {tried} times continuously, '
                      f'seems that the qsub system have some problem')
                raise e
    if 'job' not in qstat_result_string:
        # qstat print nothing, all job finished
        return []
    else:
        qstat_result_lines = qstat_result_string.split('\n')
        total_id = []
        for line in qstat_result_lines:
            # Typical 1st line of qstat
            # job-ID prior name user state submit/start at queue slots ja-task-ID
            line_id = line.strip().split(' ')[0]  # id part
            if id_set is not None:
                # only record id within id_set (submitted jobs by this project)
                if line_id in id_set:
                    total_id.append(line_id)
            else:
                # record any id if id_set is not provided
                total_id.append(line_id)
        return total_id


class _Qsubmitter:
    def __init__(self, command_file_path, working_dir, project_name, force_redo=False,
                 total_cpu=60, total_mem=500, submission_gap=SUBMISSION_GAP, qstat_gap=QSTAT_GAP):
        # prepare working_dir
        self.working_dir = working_dir
        if project_name[0].isdigit():
            raise ValueError('Project name can not start with number, qsub will fail.')
        self.project_name = project_name
        self.project_dir = self.working_dir + '/' + self.project_name + '_qsub'
        if os.path.exists(self.project_dir):
            if force_redo:
                shutil.rmtree(self.project_dir)
                os.mkdir(self.project_dir)
            else:
                print('Note: The project has been submitted before. Finished jobs will not run again.')
                print(f'If you want to redo everything, set force_redo=True')
                sys.stdout.flush()
        else:
            os.mkdir(self.project_dir)
        run(['cp', command_file_path, self.project_dir])

        # submission parameters
        self.total_cpu = total_cpu
        self.total_mem = total_mem  # GBs

        self.submission_gap = submission_gap
        self.qstat_gap = qstat_gap

        # auto detect user name
        self.user_name = run(['whoami'], stdout=PIPE, encoding='utf8').stdout.strip()
        print(f'Username: {self.user_name}')

        # prepare project dir
        self.command_file_path = command_file_path
        if not os.path.exists(self.project_dir):
            os.mkdir(self.project_dir)

        # get all Command objects
        # self.commands is a list of command objects
        self.commands = self._parse_command_file()

        # initiate running stats
        self.running_commands = []
        self.running_cpu = 0
        self.running_mem = 0  # GBs
        self.submitted_qsub_id = set()
        self.submission_fail_commands = []
        self.finished_commands = []
        self.finish_count = 0

        # add a signal file in project dir
        self.alive = True
        self.core_signal_path = self.project_dir + '/RUNNING_SIGNAL.json'
        with open(self.core_signal_path, 'w') as f:
            core_signal_dict = {'total_cpu': self.total_cpu,
                                'total_mem': self.total_mem,
                                'alive': self.alive}
            self.core_signal_dict = core_signal_dict
            json.dump(core_signal_dict, f)
        print(f"The CORE SIGNAL json file is generated: {self.core_signal_path}.\n"
              f"You can modify value (not key) in this json file "
              f"to change hot-coding parameters of the queue.\n"
              f"IMPORTANT: If 'alive' in this file is not true, qsub will stop submit "
              f"new jobs and exit until all current jobs finished.")

        # submit jobs
        self.submit()

        # get job reports
        self.get_reports()
        return

    def _parse_command_file(self):
        with open(self.command_file_path) as command:
            command_dict_list = json.load(command)
            obj_list = []
            for n, command_dict in enumerate(command_dict_list):
                obj_list.append(_Command(command_dict=command_dict,
                                         unique_id=f'{self.project_name}_{n}',
                                         working_dir=self.working_dir,
                                         project_dir=self.project_dir))
        return obj_list

    def submit(self):
        # job submission process:
        for command_obj in self.commands:
            if command_obj.finish:
                # TODO: check return code, whether to skip failed jobs
                # happens when command_obj.status_path exist
                self.finished_commands.append(command_obj)
                continue
            command_cpu = int(command_obj.qsub_parameter['pe smp'])
            command_mem = int(command_obj.qsub_parameter['l h_vmem'].strip('G'))
            future_cpu = self.running_cpu + command_cpu
            future_mem = self.running_mem + command_mem
            # check both cpu and mem
            if (future_cpu > self.total_cpu) or (future_mem > self.total_mem):
                # submitted enough jobs, wait until running_cpu decrease
                temp_gap = self.qstat_gap
                while True:
                    time.sleep(temp_gap)
                    # wait longer and longer if the job is stuck
                    if temp_gap < 300:
                        temp_gap += 5
                    self.check_running()
                    future_cpu = self.running_cpu + command_cpu
                    future_mem = self.running_mem + command_mem
                    if (future_cpu <= self.total_cpu) and (future_mem <= self.total_mem):
                        break
            command_obj.submit()
            # gather submitted id and obj
            if command_obj.submission_fail:
                self.submission_fail_commands.append(command_obj)
            else:
                self.running_commands.append(command_obj)
                self.submitted_qsub_id.add(command_obj.qsub_id)
            self.running_cpu += command_cpu
            self.running_mem += command_mem
            time.sleep(self.submission_gap)
            if not self.alive:
                print(f'Stop submitting jobs because CORE SIGNAL "alive" is not true. '
                      f'CORE SIGNAL dict path {self.core_signal_path}')
                break

        # final job check:
        temp_gap = self.qstat_gap
        while True:
            time.sleep(temp_gap)
            if temp_gap < 300:
                temp_gap += 5
            self.check_running()
            if self.running_cpu == 0:
                break  # all job finished
        return

    def check_running(self):
        # In this function, check
        # 1. running jobs, their info and resource consumption
        # 2. core_signal_dict, to adjust global control
        print('check running job: ', end='')
        cur_running_qsub_id = _get_running_job_id_qstat(user_name=self.user_name,
                                                        id_set=self.submitted_qsub_id)

        # Part 1: check every running obj
        cur_running_cpu = 0
        cur_running_mem = 0  # GBs
        cur_running_command = []
        for command_obj in self.running_commands:
            if command_obj.qsub_id not in cur_running_qsub_id:
                # print(self.submitted_qsub_id)
                # print(cur_running_qsub_id)
                # print(command_obj.qsub_id, command_obj.unique_id, 'finishing')

                # command have finished
                command_obj.finish = True
                command_obj.check_output_log()
                command_obj.write_status()
                self.finish_count += 1
                self.finished_commands.append(command_obj)
                self.submitted_qsub_id.remove(command_obj.qsub_id)
            else:
                cur_running_command.append(command_obj)
                cur_running_cpu += int(command_obj.qsub_parameter['pe smp'])
                cur_running_mem += int(command_obj.qsub_parameter['l h_vmem'].strip('G'))

        print(f'{self.finish_count} job finished in this submission.')
        sys.stdout.flush()

        self.running_commands = cur_running_command
        # update running cpu
        self.running_cpu = cur_running_cpu
        self.running_mem = cur_running_mem

        # part 2: core_signal_dict
        try:
            with open(self.core_signal_path) as f:
                self.core_signal_dict = json.load(f)
                if self.total_mem != self.core_signal_dict['total_mem']:
                    print(f"Change the queue total MEM to {self.core_signal_dict['total_mem']}")
                    self.total_mem = self.core_signal_dict['total_mem']
                if self.total_cpu != self.core_signal_dict['total_cpu']:
                    print(f"Change the queue total CPU to {self.core_signal_dict['total_cpu']}")
                    self.total_cpu = self.core_signal_dict['total_cpu']
                self.alive = self.core_signal_dict['alive']
        except OSError:
            print(f'Warning: Unable to load CORE SIGNAL json file, '
                  f'will not update core signal dict. '
                  f'This autogenerated file is a dict json, '
                  f'should exist in {self.command_file_path}.')
        return

    def get_reports(self):
        command_records = []
        for command in self.commands:
            if os.path.exists(command.status_path):
                status_series = pd.Series(command.check_submitted_status())
            else:
                status_series = pd.Series({})
            status_series['unique_id'] = command.unique_id
            status_series['submission_fail'] = command.submission_fail
            command_records.append(status_series)

        status_df = pd.DataFrame(command_records)
        status_df.set_index('unique_id', inplace=True)
        status_df.to_csv(self.project_dir + '/summary_status.tsv', sep='\t')

        print('Job submission status:')
        failed = status_df['submission_fail'].sum()
        print(f'{failed} / {status_df.shape[0]} job(s) failed in submission.')

        print('Job exit status:')
        for code, sub_df in status_df[~status_df['submission_fail']].groupby('return_code'):
            running_time = sub_df.duration_second.mean()
            time_str = str(datetime.timedelta(seconds=running_time))
            print(f'{sub_df.shape[0]} job(s) end with return code {code:.0f}, '
                  f'the average running time is {time_str}')
        return


class _Command:
    """Command object for single command script file"""

    def __init__(self, command_dict, unique_id, working_dir, project_dir):
        self.project_dir = project_dir
        self.unique_id = unique_id  # unique id is given by submitter, qsub_id is given by qsub system

        # command status
        self.submitted = False
        self.qsub_id = None
        self.submission_fail = False
        self.finish = False  # if submitted is True and not in qstat
        self.status_path = f'{self.project_dir}/{unique_id}.json'
        self.script_path = f'{self.project_dir}/{unique_id}.sh'
        self.error_path = f'{self.project_dir}/{unique_id}.error.log'
        self.output_path = f'{self.project_dir}/{unique_id}.output.log'
        self.start_time = None  # will be datetime.datetime after check_output_log
        self.end_time = None  # will be datetime.datetime after check_output_log
        self.duration_second = .0
        self.return_code = -99

        # prepare command dict
        self.command_dict = _default_command_dict(name=self.unique_id,
                                                  error_path=self.error_path,
                                                  output_path=self.output_path,
                                                  working_dir=working_dir)
        self.command_dict.update(**command_dict)
        for k, v in self.command_dict.items():
            if v is None:
                raise ValueError(f'{k} is None in command dict')

        # separate command and qsub parameter
        self.command = self.command_dict.pop('command')
        self.qsub_parameter = self.command_dict
        self.command_dict = None

        # make command sh file
        self._make_script_file()
        # if self.submit_status exist, self.submitted, self.finish will be true
        self.check_submitted_status()
        return

    def _make_script_file(self):
        """Make command file, add some echos to print start, end and return code"""
        if os.path.exists(self.script_path):
            return
        with open(self.script_path, 'w') as sh:
            sh.write("#!/bin/bash\n")
            for k, v in self.qsub_parameter.items():
                if k.startswith('l '):
                    sh.write(f'#$ -{k}={v}\n')
                else:
                    sh.write(f'#$ -{k} {v}\n')
            sh.write(f'echo JOB_CMD_START {self.unique_id} $(date +"%H:%M:%S-%D")\n')
            sh.write(self.command + '\n')
            sh.write(f'echo JOB_CMD_RETURN_CODE {self.unique_id} $?\n')  # check command return code
            sh.write(f'echo JOB_CMD_END {self.unique_id} $(date +"%H:%M:%S-%D")\n')
        return

    def submit(self):
        # there is no force submit option, because the force_redo is under control by whole project
        if os.path.exists(self.status_path):
            # if status file exist, just read the previous status, not submit job
            self.check_submitted_status()
            return
        else:
            job_id_pattern = re.compile(r'(?<=Your job )\d+')
            return_obj = run(['qsub', self.script_path], stdout=PIPE, encoding='utf8')
            print(f'Submitted {self.unique_id}, stdout: {return_obj.stdout.strip()}')
            try:
                self.qsub_id = job_id_pattern.search(return_obj.stdout).group()
            except AttributeError:
                self.submission_fail = True
            self.submitted = True
            return

    def check_submitted_status(self):
        if not os.path.exists(self.status_path):
            # do noting
            return None

        with open(self.status_path) as f:
            status_dict = json.load(f)

        self.submitted = True
        self.finish = True
        self.qsub_id = status_dict['qsub_id']
        self.start_time = datetime.datetime.strptime(status_dict['start_time'],
                                                     "%H:%M:%S-%m/%d/%y")
        self.end_time = datetime.datetime.strptime(status_dict['end_time'],
                                                   "%H:%M:%S-%m/%d/%y")
        self.duration_second = status_dict['duration_second']
        self.return_code = status_dict['return_code']
        return status_dict

    def write_status(self):
        # store the qsub id and job status.
        with open(self.status_path, 'w') as f:
            json.dump({'qsub_id': self.qsub_id,
                       'start_time': self.start_time.strftime(format='%H:%M:%S-%m/%d/%y'),
                       'end_time': self.end_time.strftime(format='%H:%M:%S-%m/%d/%y'),
                       'duration_second': self.duration_second,
                       'return_code': self.return_code}, f)
        return

    def check_output_log(self):
        if not self.finish:
            # job not finish, do nothing
            return
        with open(self.output_path) as f:
            for line in f:
                if (self.unique_id in line) and line.startswith('JOB_CMD'):
                    if line.startswith('JOB_CMD_START'):
                        self.start_time = datetime.datetime.strptime(line.strip().split(' ')[-1], "%H:%M:%S-%m/%d/%y")
                    elif line.startswith('JOB_CMD_END'):
                        self.end_time = datetime.datetime.strptime(line.strip().split(' ')[-1], "%H:%M:%S-%m/%d/%y")
                    elif line.startswith('JOB_CMD_RETURN_CODE'):
                        self.return_code = int(line.split(' ')[-1])
                    else:
                        continue
        if (self.return_code == -99) or (self.start_time is None) or (self.end_time is None):
            print('Job {self.qsub_id} return code: {self.return_code} (default init is -99)')
            print('Job {self.qsub_id} start time: {self.start_time} (default is None)')
            print('Job {self.qsub_id} end time: {self.end_time} (default is None)')
            raise ValueError(f'Output log for job {self.qsub_id} {self.unique_id} is not complete')
        self.duration_second = (self.end_time - self.start_time).total_seconds()
        return


def qsub(command_file_path, working_dir, project_name, wait_until=None,
         total_cpu=60, total_mem=500, force_redo=False,
         submission_gap=SUBMISSION_GAP, qstat_gap=QSTAT_GAP):
    if isinstance(command_file_path, str):
        command_file_path = [command_file_path]

    if wait_until is not None:
        if isinstance(wait_until, str):
            wait_until = {wait_until}
        elif isinstance(wait_until, list):
            wait_until = set(wait_until)
        else:
            raise TypeError(f'wait_until should either be str or list, '
                            f'got {type(wait_until)}')

        snap = 60
        while True:
            user_name = run(['whoami'], stdout=PIPE, encoding='utf8').stdout.strip()
            cur_running_qsub_id = _get_running_job_id_qstat(user_name, id_set=wait_until)
            if len(cur_running_qsub_id) == 0:
                break
            print(f'Still waiting {len(cur_running_qsub_id)} jobs to finish ', end='')
            if len(cur_running_qsub_id) < 10:
                job_ids = ' '.join(cur_running_qsub_id)
                print(f'Their qsub ids are: {job_ids}')
            time.sleep(snap)
            snap += 60
            snap = min(snap, 600)

    # for multiple command files, run one by one
    for i, command_file in enumerate(command_file_path):
        print(f'Execute command file: {command_file}')
        _project_name = project_name
        if len(command_file_path) != 1:
            _project_name += f'_{i}'
        submitter = _Qsubmitter(command_file_path=command_file,
                                working_dir=working_dir,
                                project_name=_project_name,
                                total_cpu=total_cpu,
                                total_mem=total_mem,
                                submission_gap=submission_gap,
                                qstat_gap=qstat_gap,
                                force_redo=force_redo)
        print(f'{len(submitter.commands)} jobs finished.')
    return
