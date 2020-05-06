"""
Package file structure:

test/
    pytest data and code
mapping/
    dealing with snmc-seq2 hdf5 mapping, from fastq to allc, and generate QC stats.
deprecated_plot/
    all plotting functions related to any step of analysis
tools/
    function for analysis and manipulating text based files (allc, bed, bigwig etc.).
local/
    unstable code that only used for local analysis, may be changed or removed in any time.
__main__.py
    enter point of yap CLI
__init__.py
    enter point of python API

# TODO add api for useful functions, like scanpy
"""

__version__ = '0.1.2'
