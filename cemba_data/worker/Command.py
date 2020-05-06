import subprocess
import shlex

COMMAND_STATUS = [
    'PD',  # Pending
    'R',  # Running
    'S',  # Success
    'E',  # Error
    'SKIP'  # Deliberately skipped due to dependent command is skipped or failed
]


class Command:
    """
    Command is a wrapper of subprocess.run, in addition to subprocess.run,
    it check whether the dependent commands finished, and manage status of this command

    Lifespan of a command obj
    1. init with command and dependency
    2. pre_run check
    3. run
    4. post_run check
    5. save to some archive
    """

    def __init__(self, command, dependent_commands=None, core=1, mem=3):
        self.command = command
        self.core = core
        self.mem = mem
        self._status = 'PD'
        self._runnable = False

        if dependent_commands is None:
            self._dependent_list = []
        else:
            self._dependent_list = dependent_commands

        self.stdout = ''
        self.stderr = ''
        self.return_code = ''
        self.stats = {}  # save other stats file and information, manage this by post_run
        return

    @property
    def status(self):
        return self._status

    @status.setter
    def status(self, value):
        value = value.upper()
        if value not in COMMAND_STATUS:
            raise ValueError(f'Status of a command can only be {COMMAND_STATUS}, got {value}')
        self._status = value
        return

    def pre_run(self):
        """Rewrite this for specific validations pre run,
        this function may have more complex checks and change the command.status"""
        return

    def post_run(self):
        """Rewrite this for specific validations post run,
        this function may have more post checks, clean up temp file or collect status file
        and change the command.status"""
        return

    def runnable(self):
        if not self._runnable:
            if len(self._dependent_list) != 0:
                status_set = set([c.status for c in self._dependent_list])
                if ('PD' in status_set) or ('R' in status_set):
                    # this command still need to wait for dependent commands
                    self.status = 'PD'
                    self._runnable = False
                elif ('E' in status_set) or ('SKIP' in status_set):
                    self.status = 'SKIP'
                    self._runnable = True
                else:
                    # only "S" means all dependent commands successfully executed
                    self.status = 'R'
                    self._runnable = True
            else:
                # no dependency
                self.status = 'R'
                self._runnable = True
        return self._runnable

    def run(self):
        """
        Check dependency and run, return status
        """
        self.runnable()
        if not self._runnable:
            raise ValueError(f'The command is not able to run yet, because the dependent commands have not finish.')

        # check status before run, this can happen when self.pre_run have some issue
        if self.status in ['E', 'SKIP']:
            return self.status

        # run the command
        try:
            obj = subprocess.run(shlex.split(self.command),
                                 stderr=subprocess.PIPE,
                                 stdout=subprocess.PIPE,
                                 encoding='utf8',
                                 check=True)
        except subprocess.CalledProcessError as obj:
            self.status = 'E'
        else:
            self.status = 'S'

        self.return_code = obj.returncode
        self.stderr = obj.stderr
        self.stdout = obj.stdout

        return self.status
