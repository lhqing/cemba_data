"""
- read and write allc file
- parallel writing
- read region from allc
- separate process gives better performance


This file is modified from xopen 0.3.4
Here is the licence

Copyright (c) 2010-2018 Marcel Martin <mail@marcelm.net>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

from __future__ import print_function, division, absolute_import

import gzip
import sys
import os
import time
import codecs
from subprocess import Popen, PIPE, run

try:
    run(['pigz', '--version'],
        stdout=PIPE,
        stderr=PIPE)
    PIGZ = True
except OSError:
    PIGZ = False

try:
    run(['bgzip', '--version'],
        stdout=PIPE,
        stderr=PIPE)
    BGZIP = True
except OSError:
    BGZIP = False


class Closing(object):
    """
    Inherit from this class and implement a close() method to offer context
    manager functionality.
    """

    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __del__(self):
        try:
            self.close()
        except OSError:
            pass


class PipedGzipWriter(Closing):
    """
    Write gzip-compressed files by running an external gzip or pigz process and
    piping into it. On Python 2, this is faster than using gzip.open(). On
    Python 3, it allows to run the compression in a separate process and can
    therefore also be faster.
    """

    def __init__(self, path, mode='wt', compresslevel=6, threads=1):
        """
        mode -- one of 'w', 'wt', 'wb', 'a', 'at', 'ab'
        compresslevel -- gzip compression level
        threads (int) -- number of pigz threads (None means to let pigz decide)
        """
        if mode not in ('w', 'wt', 'wb', 'a', 'at', 'ab'):
            raise ValueError("Mode is '{0}', but it must be 'w', 'wt', 'wb', 'a', 'at' or 'ab'".format(mode))

        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, mode)
        self.closed = False
        self.name = path

        kwargs = dict(stdin=PIPE, stdout=self.outfile, stderr=self.devnull)
        # Setting close_fds to True in the Popen arguments is necessary due to
        # <http://bugs.python.org/issue12786>.
        # However, close_fds is not supported on Windows. See
        # <https://github.com/marcelm/cutadapt/issues/315>.
        if sys.platform != 'win32':
            kwargs['close_fds'] = True

        # only apply to gzip and pigz, bgzip use -l
        if 'w' in mode and compresslevel != 6:
            extra_args = ['-' + str(compresslevel)]
        else:
            extra_args = []

        if 'b' not in mode:
            encoding = None
        else:
            encoding = 'utf8'
        kwargs['encoding'] = encoding

        try:
            if BGZIP:
                bgzip_args = ['bgzip']
                if threads is not None and threads > 0:
                    bgzip_args += ['-@', str(threads), '-l', str(compresslevel)]
                self.process = Popen(bgzip_args, **kwargs)
                self.program = 'bgzip'
            elif PIGZ:
                pigz_args = ['pigz']
                if threads is not None and threads > 0:
                    pigz_args += ['-p', str(threads)]
                self.process = Popen(pigz_args + extra_args, **kwargs)
                self.program = 'pigz'
            else:
                # pigz not found, try regular gzip
                self.process = Popen(['gzip'] + extra_args, **kwargs)
                self.program = 'gzip'
        except OSError:
            self.outfile.close()
            self.devnull.close()
            raise
        self._file = codecs.getwriter('utf-8')(self.process.stdin)

    def write(self, arg):
        self._file.write(arg)

    def close(self):
        self.closed = True
        self._file.close()
        return_code = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if return_code != 0:
            raise IOError(f"Output {self.program} process terminated with exit code {return_code}")


class PipedGzipReader(Closing):
    # decompression can't be parallel even in pigz, so there is not thread/cpu parameter
    def __init__(self, path, region=None, mode='r'):
        if mode not in ('r', 'rt', 'rb'):
            raise ValueError(f"Mode is {mode}, but it must be 'r', 'rt' or 'rb'")
        if 'b' not in mode:
            encoding = 'utf8'
        else:
            encoding = None

        if region is None:
            self.process = Popen(['gzip', '-cd', path],
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)
        else:
            self.process = Popen(['tabix', path, region],
                                 stdout=PIPE,
                                 stderr=PIPE,
                                 encoding=encoding)

        self.name = path
        self._file = self.process.stdout
        self._stderr = self.process.stderr
        self.closed = False
        # Give gzip a little bit of time to report any errors
        # (such as a non-existing file)
        time.sleep(0.01)
        self._raise_if_error()

    def close(self):
        self.closed = True
        return_code = self.process.poll()
        if return_code is None:
            # still running
            self.process.terminate()
        self._raise_if_error()

    def __iter__(self):
        for line in self._file:
            yield line
        self.process.wait()
        self._raise_if_error()

    def readline(self):
        return self._file.readline()

    def _raise_if_error(self):
        """
        Raise IOError if process is not running anymore and the
        exit code is nonzero.
        """
        return_code = self.process.poll()
        if return_code is not None and return_code != 0:
            message = self._stderr.read().strip()
            raise IOError(message)

    def read(self, *args):
        data = self._file.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data


def open_gz(filename, mode='r', compresslevel=6, threads=1, region=None):
    if 'r' in mode:
        try:
            return PipedGzipReader(filename, region=region, mode=mode)
        except OSError:
            # gzip not installed
            return gzip.open(filename, mode)
    else:
        try:
            return PipedGzipWriter(filename, mode, compresslevel, threads=threads)
        except OSError:
            return gzip.open(filename, mode, compresslevel=compresslevel)


def open_allc(filename, mode='r', compresslevel=3, threads=1,
              region=None):
    """
    A replacement for the "open" function that can also open files that have
    been compressed with gzip, bzip2 or xz. If the filename is '-', standard
    output (mode 'w') or input (mode 'r') is returned.

    The file type is determined based on the filename: .gz is gzip, .bz2 is bzip2 and .xz is
    xz/lzma.

    When writing a gzip-compressed file, the following methods are tried in order to get the
    best speed 1) using a pigz (parallel gzip) subprocess; 2) using a gzip subprocess;
    3) gzip.open. A single gzip subprocess can be faster than gzip.open because it runs in a
    separate process.

    Uncompressed files are opened with the regular open().

    mode can be: 'rt', 'rb', 'at', 'ab', 'wt', or 'wb'. Also, the 't' can be omitted,
    so instead of 'rt', 'wt' and 'at', the abbreviations 'r', 'w' and 'a' can be used.

    threads is the number of threads for pigz. If None, then the pigz default is used.
    multi-thread only apply to writer, reader (decompression) can't be paralleled
    """
    if mode in ('r', 'w', 'a'):
        mode += 't'
    if mode not in ('rt', 'rb', 'wt', 'wb', 'at', 'ab'):
        raise ValueError("mode '{0}' not supported".format(mode))
    if not isinstance(filename, str):
        raise ValueError("the filename must be a string")
    if compresslevel not in range(1, 10):
        raise ValueError("compresslevel must be between 1 and 9")
    if region is not None:
        # unzipped file
        if not filename.endswith('gz'):
            raise ValueError('File must be compressed by bgzip to use region query.')
        # normal gzipped file
        if not has_tabix(filename):
            raise ValueError(f'Tried inspect {filename}, '
                             'File is compressed by normal gzip, but region query only apply to bgzip')

        if not os.path.exists(filename + '.tbi'):
            raise FileNotFoundError('region query provided, but .tbi index not found')

    if filename.endswith('gz'):
        return open_gz(filename, mode, compresslevel, threads, region=region)
    else:
        return open(filename, mode)


def rezip_use_bgzip(filename):
    if is_bgzip(filename):
        return

    tmp_filename = filename + '_rezipping'
    if os.path.exists(tmp_filename):
        raise FileExistsError(f'Temp file {tmp_filename} exist')
    run(['mv', filename, tmp_filename], check=True)
    reader = Popen(['gzip', '-cd', tmp_filename],
                   stdout=PIPE).stdout
    try:
        with open(filename, 'w') as f:
            writer = Popen(['bgzip'], stdin=reader, stdout=f)
            writer.wait()
        run(['rm', '-f', tmp_filename])
    except OSError as e:
        run(['mv', tmp_filename, filename], check=True)
        raise e
    return


def is_bgzip(filename):
    # TODO This function is not perfect, figure out a more substantial way to determine bgzip
    # a bad example is here, its bgziped, but file provide something strange
    # /gale/raidix/rdx-4/CEMBA_RS1/4C/CEMBA180417_4C/allc/
    # allc_180605_CEMBA_mm_P56_P63_4C_CEMBA180417_4C_7_CEMBA180417_4C_8_C8_AD001_indexed.tsv.gz
    print('Warning: is_bgzip(), This function is not perfect, figure out a more substantial way to determine bgzip')
    if not filename.endswith('gz'):
        return False

    file_type_descript = run(['file', filename],
                             stdout=PIPE, encoding='utf8').stdout
    # for bgzip format, I don't have a better idea to check
    # this "gzip compressed data, extra field" description is different from normal gzip
    # for normal gzip, result looks like:
    # test.tsv.gz: gzip compressed data, was "test.tsv", from Unix, last modified: Fri Jan  4 17:17:31 2019
    if 'gzip compressed data, extra field' in file_type_descript:
        return True
    else:
        return False


def has_tabix(filename):
    tabix_path = filename + '.tbi'
    if os.path.exists(tabix_path):
        return True
    else:
        return False
