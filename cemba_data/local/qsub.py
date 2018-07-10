from subprocess import run, PIPE
import os
import json
import configparser
import datetime
import re
import time

qsub_config = configparser.ConfigParser()
qsub_config.read(os.path.dirname(__file__) + '/config_qsub.ini')


def default_command_dict(name, error_path, output_path):
    command_dict = {'command': None,
                    '-N': name,
                    '-V': '',
                    '-pe smp': qsub_config['QSUB_DEFAULT']['CPU'],
                    '-l h_rt': qsub_config['QSUB_DEFAULT']['H_RT'],
                    '-l s_rt': qsub_config['QSUB_DEFAULT']['S_RT'],
                    '-wd': qsub_config['QSUB_DEFAULT']['WD'],
                    '-e': error_path,
                    '-o': output_path}
    return command_dict


def check_qstat(id_set=None):
    qstat_result_string = run(['qstat', '-u', qsub_config['USER']['USER_NAME']],
                              stdout=PIPE, encoding='utf8').stdout
    if 'job-ID' not in qstat_result_string:
        return []  # qstat print nothing, all job finished
    else:
        qstat_result_lines = qstat_result_string.split('\n')
        total_id = []
        for line in qstat_result_lines:
            line_id = line.split(' ')[0]
            if id_set is not None:
                if line_id in id_set:
                    total_id.append(line_id)
            else:
                total_id.append(line_id)
        return total_id


class Qsubmitter:
    def __init__(self, command_file_path, project_name,
                 auto_submit=False, qacct=False, resubmit=None,
                 total_cpu=None, submission_gap=None, qstat_gap=None):
        self.command_file_path = command_file_path
        self.project_name = project_name
        self.project_dir = qsub_config['QSUB_DEFAULT']['JOB_DIR'] + '/' + project_name
        if not os.path.exists(self.project_dir):
            os.mkdir(self.project_dir)

        self.commands = self.parse_command_file()
        self.running_commands = []
        self.running_cpu = 0
        self.submitted_qsub_id = set()
        self.submission_fail_commands = []
        self.finished_commands = []  # modify this with check_command_finish_status
        self.qsub_failed_commands = []  # modify this with check_command_finish_status
        self.cmd_failed_commands = []  # modify this with check_command_finish_status
        self.check_qacct = qacct  # check qacct is slow if the log is too large, only check it in debugging
        self.finish_count = 0
        if total_cpu is None:
            self.total_cpu = int(qsub_config['RUNNING_DEFAULT']['TOTAL_CPU'])
        else:
            self.total_cpu = total_cpu
        if submission_gap is None:
            self.submission_gap = int(qsub_config['RUNNING_DEFAULT']['SUBMISSION_GAP'])
        else:
            self.submission_gap = submission_gap
        if qstat_gap is None:
            self.qstat_gap = int(qsub_config['RUNNING_DEFAULT']['QSTST_GAP'])
        else:
            self.qstat_gap = qstat_gap

        if auto_submit:
            self.submit(resubmit)
        return

    def parse_command_file(self):
        with open(self.command_file_path) as command:
            command_dict_list = json.load(command)
            obj_list = []
            n = 0
            for command_dict in command_dict_list:
                obj_list.append(Command(command_dict=command_dict,
                                        unique_id=f'{self.project_name}_{n}',
                                        project_dir=self.project_dir))
                n += 1
        return obj_list

    def submit(self, resubmit=None):
        # job submission process:
        for command_obj in self.commands:
            # judge if command has been submitted
            if command_obj.submitted:
                if resubmit is None:
                    # print('Skip command:', command_obj.unique_id)
                    self.check_command_finish_status(command_obj)
                    continue  # not resubmit command
                elif resubmit == 'fail':
                    # only resubmit failed command
                    if not command_obj.qsub_failed and not command_obj.cmd_failed:
                        # both a success, the job is success
                        # print('Skip command:', command_obj.unique_id)
                        self.check_command_finish_status(command_obj)
                        continue
                    else:
                        print('Resubmit failed command:', command_obj.unique_id)
                elif resubmit == 'all':
                    # resubmit all command
                    pass
                else:
                    raise ValueError('resubmit only allow None, fail or all. Got', resubmit)

            command_cpu = int(command_obj.qsub_parameter['-pe smp'])
            if self.running_cpu + command_cpu > self.total_cpu:
                while True:
                    time.sleep(self.qstat_gap)
                    self.check_running_cpu()
                    if self.running_cpu <= self.total_cpu:
                        break
            # print('Submit job:', command_obj.unique_id)
            command_obj.submit()
            # gather submitted id and obj
            if command_obj.submission_fail:
                self.submission_fail_commands.append(command_obj)
            else:
                self.running_commands.append(command_obj)
                self.submitted_qsub_id.add(command_obj.qsub_id)
            self.running_cpu += command_cpu
            time.sleep(self.submission_gap)

        # final job check:
        while True:
            time.sleep(self.qstat_gap)
            self.check_running_cpu()
            if self.running_cpu == 0:
                break  # all job finished
        self.get_reports()
        return

    def check_running_cpu(self):
        # TODO: check result in the .error.log to see if there is any error
        print('check running job')
        cur_running_qsub_id = check_qstat(self.submitted_qsub_id)
        # check every running obj
        cur_running_cpu = 0
        cur_running_command = []
        for command_obj in self.running_commands:
            if command_obj.qsub_id not in cur_running_qsub_id:
                # command have finished
                self.finish_count += 1

                if self.check_qacct:
                    print(command_obj.unique_id, 'qacct')
                    if command_obj.query_qacct() != 0:
                        raise ValueError('Command not finished but disappeared in qstat')
                else:
                    command_obj.not_query_qacct()
                self.check_command_finish_status(command_obj)
            else:
                cur_running_command.append(command_obj)
                cur_running_cpu += int(command_obj.qsub_parameter['-pe smp'])
        print(f'{self.finish_count} / {len(self.commands)} job finished in this submission.')
        self.running_commands = cur_running_command
        # update running cpu
        self.running_cpu = cur_running_cpu
        return

    def check_command_finish_status(self, command_obj):
        self.finished_commands.append(command_obj)
        if self.check_qacct:
            if command_obj.qsub_failed:
                self.qsub_failed_commands.append(command_obj)
            if command_obj.cmd_failed:
                self.cmd_failed_commands.append(command_obj)
        return

    def get_reports(self):
        cur_time = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        with open(self.project_dir + f'/submit_summary.{cur_time}.txt', 'w') as summary:
            summary.write(f'Total Command: {len(self.commands)}\n')
            summary.write(f'Submission Failed Command: {len(self.submission_fail_commands)}\n')

            if self.check_qacct:
                summary.write(f'Qsub Failed Command: {len(self.qsub_failed_commands)}\n')
                if len(self.qsub_failed_commands) > 0:
                    with open(self.project_dir + '/qsub_failed_commands.txt', 'w') as f:
                        f.write('\n'.join(self.qsub_failed_commands))
                    summary.write(f'See qsub failed command id in {self.project_dir}/qsub_failed_commands.txt\n')

                summary.write(f'Command Failed Command: {len(self.cmd_failed_commands)}\n')
                if len(self.cmd_failed_commands) > 0:
                    with open(self.project_dir + '/cmd_failed_commands.txt', 'w') as f:
                        f.write('\n'.join(self.cmd_failed_commands))
                    summary.write(f'See qsub failed command id in {self.project_dir}/cmd_failed_commands.txt\n')
            else:
                summary.write('check_qacct = False, '
                              'if want to get the finish code for each job, '
                              'set qacct=True in Qsubmitter')
        return


class Command:
    def __init__(self, command_dict, unique_id, project_dir,
                 error_path=None, output_path=None):
        self.project_dir = project_dir
        self.unique_id = unique_id

        # command status
        self.submitted = False
        self.qsub_id = None
        self.submission_fail = False
        self.finish = False  # if submitted is True and not in qstat
        self.qsub_failed = None  # check from qacct_status['failed'], fail due to qsub
        self.cmd_failed = None  # check from qacct_status['exit_status'], fail due to cmd
        self.qacct_status = None  # dict of qacct results
        self.qacct_status_path = f'{self.project_dir}/{unique_id}.qacct.json'
        self.script_path = f'{self.project_dir}/{unique_id}.sh'
        if error_path is None:
            self.error_path = f'{self.project_dir}/{unique_id}.error.log'
        if output_path is None:
            self.output_path = f'{self.project_dir}/{unique_id}.output.log'

        # prepare command dict
        self.command_dict = default_command_dict(name=self.unique_id,
                                                 error_path=self.error_path,
                                                 output_path=self.output_path)
        self.command_dict.update(**command_dict)
        for k, v in self.command_dict.items():
            if v is None:
                raise ValueError(f'{k} is None in command dict')

        # separate command and qsub parameter
        self.command = self.command_dict.pop('command')
        self.qsub_parameter = self.command_dict
        self.command_dict = None

        # make command sh file
        self.make_script_file()

        # modify command status if it has been submitted and finished
        if os.path.exists(self.qacct_status_path):
            with open(self.qacct_status_path) as f:
                status_dict = json.load(f)
                if len(status_dict) > 1:
                    self.qacct_status = status_dict
                    self.qsub_id = self.qacct_status['jobnumber']
                    self.judge_qacct_status()
                else:  # didn't do qacct
                    self.qsub_id = status_dict['qsub_id']
                    # self.qacct_status is still None,
                    # should run query_qacct to get it
            self.submitted = True
            self.finish = True
        return

    def make_script_file(self):
        with open(self.script_path, 'w') as sh:
            sh.write("#!/bin/bash\n")
            for k, v in self.qsub_parameter.items():
                if k[:2] == '-l':
                    sh.write(f'#$ {k}={v}\n')
                else:
                    sh.write(f'#$ {k} {v}\n')
            sh.write(self.command)
        return

    def submit(self):
        job_id_pattern = re.compile(r'(?<=Your job )\d+')
        return_obj = run(['qsub', self.script_path], stdout=PIPE, encoding='utf8')
        try:
            self.qsub_id = job_id_pattern.search(return_obj.stdout).group()
        except AttributeError:
            self.submission_fail = True
        self.submitted = True

        # remove previous qacct status
        self.finish = False
        self.qsub_failed = None
        self.cmd_failed = None
        self.qacct_status = None
        return

    def query_qacct(self):
        if self.qacct_status is not None:
            return 0
        time.sleep(10)
        qstat_result = run(['qstat', '-j', self.qsub_id], stdout=PIPE, stderr=PIPE, encoding='utf8')
        if 'Following jobs do not exist' in qstat_result.stderr and self.submitted:
            self.finish = True
        else:
            return 1  # not finish
        qacct_result = run(['qacct', '-j', self.qsub_id, '-o', qsub_config['USER']['USER_NAME']],
                           stdout=PIPE, encoding='utf8')
        status_dict = {}
        qacct_result_lines = qacct_result.stdout.split('\n')[1:-1]
        for line in qacct_result_lines:
            k = line[:13].strip()  # how to separate lines without \t???
            v = line[13:].strip()
            status_dict[k] = v
        self.qacct_status = status_dict
        self.judge_qacct_status()
        with open(self.qacct_status_path, 'w') as f:
            json.dump(self.qacct_status, f)
        return 0

    def not_query_qacct(self):
        self.finish = True
        # skip qacct, only store the qsub id.
        with open(self.qacct_status_path, 'w') as f:
            json.dump({'qsub_id': self.qsub_id}, f)

    def judge_qacct_status(self):
        if self.qacct_status['failed'] == '0':
            self.qsub_failed = False
        else:
            self.qsub_failed = True
        if self.qacct_status['exit_status'] == '0':
            self.cmd_failed = False
        else:
            self.cmd_failed = True
        return
