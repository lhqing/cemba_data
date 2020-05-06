"""
The local worker execute a list of commands locally
It needs to know:
1. total core and MEM to use
2. commands information

It does:
1. check command dependency from a center DB
2. execute commands in parallel under total limit
3. check and save command return status to a center DB
"""
from concurrent.futures import ProcessPoolExecutor, as_completed


class Executor:
    """
    Executor take a list of command objects and run the command
    """

    def __init__(self, commands, total_cpu, total_mem):
        self.pending_commands = commands
        self.total_cpu = total_cpu
        self.total_mem = total_mem

        self.runnable_command = []  # when command.runnable() return true
        self.done_command = []  # when command.status is S, E, or SKIP
        self.future_dict = {}
        return

    def refresh_runnable(self):
        _temp_pending = []
        for command in self.pending_commands:
            if command.runnable():
                self.runnable_command.append(command)
            else:
                _temp_pending.append(command)
        self.pending_commands = _temp_pending
        return

    def start(self):
        executor = ProcessPoolExecutor(self.total_cpu)
        return executor


"""
In worker process 
whenever worker want to take a task, it first check active worker dict to see whether active worker reach maximum limit,
if not, it will register itself to active worker dict and start the task.
When task finish, it will save the task result to result queue, 
remove itself from active_worker_dict and request task again
if worker got poison pill, it terminate itself

In main process
start workers, task queue, result queue, and active_worker_dict
refresh runnable, put runnable into task queue
wait for a while, refresh runnable, put runnable into task queue again, until all runnable put in
put poison pills into task queue
join queue and wait for finish
Check result queue and process final result
"""
