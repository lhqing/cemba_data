"""
SGE qsub system auto-submitter

TODO: add qacct check and better report run time and success/failure
TODO: Standardize parameter support
"""

from subprocess import run, PIPE
import os
import json
import datetime
import re
import sys
import time
import argparse


def default_command_dict(name, error_path, output_path, working_dir,
                         cpu=1, h_rt='99:99:99', s_rt='99:99:99'):
    command_dict = {'command': None,
                    '-N': name,
                    '-V': '',
                    '-pe smp': cpu,
                    '-l h_rt': h_rt,
                    '-l s_rt': s_rt,
                    '-wd': working_dir,
                    '-e': error_path,
                    '-o': output_path}
    return command_dict


def get_running_job_id_qstat(user_name, id_set):
    qstat_result_string = run(['qstat', '-u', user_name],
                              stdout=PIPE, encoding='utf8').stdout
    if 'job-ID' not in qstat_result_string:
        # print('qstat nothing')
        return []  # qstat print nothing, all job finished
    else:
        qstat_result_lines = qstat_result_string.split('\n')
        total_id = []
        for line in qstat_result_lines:
            line_id = line.strip().split(' ')[0]
            if id_set is not None:
                if line_id in id_set:
                    total_id.append(line_id)
            else:
                total_id.append(line_id)
        return total_id


class Qsubmitter:
    def __init__(self, command_file_path, working_dir, project_name,
                 total_cpu=60, submission_gap=2, qstat_gap=30):
        # prepare working_dir
        self.working_dir = working_dir
        self.project_name = project_name
        self.project_dir = self.working_dir + '/' + self.project_name + '_qsub'
        try:
            os.mkdir(self.project_dir)
        except OSError:
            print('Note: The job has been submitted before.')
            print(f'If want to resubmit everything, delete or rename this directory: {self.project_dir}')
            sys.stdout.flush()
        run(['cp', command_file_path, self.project_dir])

        # submission parameters
        self.total_cpu = total_cpu
        self.submission_gap = submission_gap
        self.qstat_gap = qstat_gap

        # auto detect user name
        self.user_name = run(['whoami'], stdout=PIPE, encoding='utf8').stdout.strip()

        # prepare project dir
        self.command_file_path = command_file_path
        if not os.path.exists(self.project_dir):
            os.mkdir(self.project_dir)

        # get all Command objects
        self.commands = self.parse_command_file()

        # initiate running stats
        self.running_commands = []
        self.running_cpu = 0
        self.submitted_qsub_id = set()
        self.submission_fail_commands = []
        self.finished_commands = []
        self.finish_count = 0

        # submit jobs
        self.submit()
        return

    def parse_command_file(self):
        with open(self.command_file_path) as command:
            command_dict_list = json.load(command)
            obj_list = []
            for n, command_dict in enumerate(command_dict_list):
                obj_list.append(Command(command_dict=command_dict,
                                        unique_id=f'{self.project_name}_{n}',
                                        working_dir=self.working_dir,
                                        project_dir=self.project_dir))
        return obj_list

    def submit(self):
        # job submission process:
        for command_obj in self.commands:
            command_cpu = int(command_obj.qsub_parameter['-pe smp'])
            if (self.running_cpu + command_cpu) > self.total_cpu:
                while True:
                    time.sleep(self.qstat_gap)
                    self.check_running_cpu()
                    if self.running_cpu <= self.total_cpu:
                        break
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
        print('check running job: ', end='')
        cur_running_qsub_id = get_running_job_id_qstat(user_name=self.user_name, id_set=self.submitted_qsub_id)
        # check every running obj
        cur_running_cpu = 0
        cur_running_command = []
        for command_obj in self.running_commands:
            if command_obj.qsub_id not in cur_running_qsub_id:
                # command have finished
                command_obj.finish = True

                # store the qsub id.
                with open(command_obj.status_path, 'w') as f:
                    json.dump({'qsub_id': command_obj.qsub_id}, f)
                self.finish_count += 1
                self.finished_commands.append(command_obj)
            else:
                cur_running_command.append(command_obj)
                cur_running_cpu += int(command_obj.qsub_parameter['-pe smp'])
        print(f'{self.finish_count} / {len(self.commands)} job finished in this submission.')
        sys.stdout.flush()
        self.running_commands = cur_running_command
        # update running cpu
        self.running_cpu = cur_running_cpu
        return

    def get_reports(self):
        cur_time = datetime.datetime.now().strftime("[%Y-%m-%d]%H:%M:%S")
        with open(self.project_dir + f'/submit_summary.{cur_time}.txt', 'w') as summary:
            summary.write(f'Total Command: {len(self.commands)}\n')
            summary.write(f'Submission Failed Command: {len(self.submission_fail_commands)}\n')
        return


class Command:
    def __init__(self, command_dict, unique_id, working_dir, project_dir):
        self.project_dir = project_dir
        self.unique_id = unique_id

        # command status
        self.submitted = False
        self.qsub_id = None
        self.submission_fail = False
        self.finish = False  # if submitted is True and not in qstat
        self.status_path = f'{self.project_dir}/{unique_id}.json'
        self.script_path = f'{self.project_dir}/{unique_id}.sh'
        self.error_path = f'{self.project_dir}/{unique_id}.error.log'
        self.output_path = f'{self.project_dir}/{unique_id}.output.log'

        # prepare command dict
        self.command_dict = default_command_dict(name=self.unique_id,
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
        self.make_script_file()
        return

    def make_script_file(self):
        with open(self.script_path, 'w') as sh:
            sh.write("#!/bin/bash\n")
            # TODO: include all parameters and set same format
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
        return


def qsub(command_file_path, working_dir, project_name, total_cpu=60, submission_gap=2, qstat_gap=30):
    submitter = Qsubmitter(command_file_path=command_file_path,
                           working_dir=working_dir,
                           project_name=project_name,
                           total_cpu=total_cpu,
                           submission_gap=submission_gap,
                           qstat_gap=qstat_gap)
    print(f'{len(submitter.commands)} jobs finished.')
    return


def qsub_register_subparser(subparser):
    parser = subparser.add_parser('qsub',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  help="General qsub helper, need to prepare a command dict file")
    parser.set_defaults(func=qsub)

    parser_req = parser.add_argument_group("Required inputs")
    parser_opt = parser.add_argument_group("Optional inputs")

    parser_req.add_argument(
        "--working_dir",
        type=str,
        required=True,
        help="Working directory of the work project"
    )

    parser_req.add_argument(
        "--project_name",
        type=str,
        required=True,
        help="Name of the work project"
    )

    parser_req.add_argument(
        "--command_file_path",
        type=str,
        required=True,
        help="Path of the command dict file"
    )

    parser_opt.add_argument(
        "--total_cpu",
        type=int,
        required=False,
        help="Total CPU in qsub list"
    )

    parser_opt.add_argument(
        "--submission_gap",
        type=int,
        required=False,
        help="Submission Gap in qsub list"
    )

    parser_opt.add_argument(
        "--qstat_gap",
        type=int,
        required=False,
        help="Qstat check gap in qsub list"
    )
    return
