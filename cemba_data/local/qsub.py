"""
SGE qsub system auto-submitter
"""

from subprocess import run, PIPE
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
    qstat_result_string = run(['qstat', '-u', user_name],
                              stdout=PIPE, encoding='utf8').stdout
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
        # TODO add running mem check too
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

        # submit jobs
        self.submit()
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
                    temp_gap += 5
                    self.check_running()
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

        # final job check:
        temp_gap = self.qstat_gap
        while True:
            time.sleep(temp_gap)
            temp_gap += 5
            self.check_running()
            if self.running_cpu == 0:
                break  # all job finished
        self.get_reports()
        return

    def check_running(self):
        print('check running job: ', end='')
        cur_running_qsub_id = _get_running_job_id_qstat(user_name=self.user_name, id_set=self.submitted_qsub_id)
        # check every running obj
        cur_running_cpu = 0
        cur_running_mem = 0  # GBs
        cur_running_command = []
        for command_obj in self.running_commands:
            if command_obj.qsub_id not in cur_running_qsub_id:
                # command have finished
                command_obj.finish = True
                command_obj.check_output_log()
                command_obj.write_status()
                self.finish_count += 1
                self.finished_commands.append(command_obj)
            else:
                cur_running_command.append(command_obj)
                cur_running_cpu += int(command_obj.qsub_parameter['pe smp'])
                cur_running_mem += int(command_obj.qsub_parameter['l h_vmem'].strip('G'))

        print(f'{self.finish_count} / {len(self.commands)} job finished in this submission.')
        sys.stdout.flush()

        self.running_commands = cur_running_command
        # update running cpu
        self.running_cpu = cur_running_cpu
        self.running_mem = cur_running_mem
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
        self.return_code = -1

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
            try:
                self.qsub_id = job_id_pattern.search(return_obj.stdout).group()
            except AttributeError:
                self.submission_fail = True
            self.submitted = True
            return

    def check_submitted_status(self):
        self.submitted = True
        with open(self.status_path) as f:
            status_dict = json.load(f)
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
        if self.return_code == -1:
            # didn't find return code line
            self.return_code = 1
        if (self.start_time is not None) and (self.end_time is not None):
            self.duration_second = (self.end_time - self.start_time).total_seconds()
        return


def qsub(command_file_path, working_dir, project_name,
         total_cpu=60, total_mem=500, force_redo=False,
         submission_gap=SUBMISSION_GAP, qstat_gap=QSTAT_GAP):
    if isinstance(command_file_path, str):
        command_file_path = [command_file_path]
    if isinstance(total_cpu, int):
        total_cpu = [total_cpu] * len(command_file_path)
    if isinstance(total_mem, int):
        total_mem = [total_mem] * len(command_file_path)

    # for multiple command files, run one by one
    for i, (command_file, cpu, mem) in enumerate(zip(command_file_path, total_cpu, total_mem)):
        print(f'Execute command file: {command_file}')
        _project_name = project_name
        if len(command_file_path) != 1:
            _project_name += f'_{i}'
        submitter = _Qsubmitter(command_file_path=command_file,
                                working_dir=working_dir,
                                project_name=_project_name,
                                total_cpu=cpu,
                                total_mem=mem,
                                submission_gap=submission_gap,
                                qstat_gap=qstat_gap,
                                force_redo=force_redo)
        print(f'{len(submitter.commands)} jobs finished.')
    return
