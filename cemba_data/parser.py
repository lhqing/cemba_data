import argparse
import sys
import datetime
import inspect

import cemba_data
from cemba_data.tools.allc import map_to_region_register_subparser
from cemba_data.local.prepare_dataset import prepare_dataset_register_subparser


def cur_time():
    return datetime.datetime.now().strftime("%Y/%m/%d-%H:%M:%S")


def get_common_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar=""
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'register_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(cemba_data.__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # execute command
    args_vars = vars(args)
    cur_func = args_vars.pop('func')
    cur_command = args_vars.pop('command')
    print(f"{cur_time()}\tExecuting {cur_command}...")
    cur_func(**args_vars)
    print(f"{cur_time()}\t{cur_command} finished.")
    return
