import shlex
import subprocess


def testing_cmd(command, expected_return_code=0):
    try:
        p = subprocess.run(shlex.split(command),
                           stderr=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           encoding='utf8',
                           check=True)
        print(p.stdout)
    except subprocess.CalledProcessError as e:
        if e.returncode == expected_return_code:
            return
        print(e.stderr)
        raise e
    return


COMMAND_TO_TEST = [
    'cutadapt --version',
    'bismark -version',
    'bowtie2 --version',
    'samtools --version',
    'tabix --version',
    'bgzip --version',
    'bedtools --version'
]


def testing_mapping_installation():
    for command in COMMAND_TO_TEST:
        testing_cmd(command)

    # picard always return 1...
    testing_cmd('picard MarkDuplicates --version', 1)
