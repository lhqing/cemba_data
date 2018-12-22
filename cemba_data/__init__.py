"""
Package file structure:

hdf5/
    dealing with HDF5 based datasets, such as netcdf file or anndata file.
mapping/
    dealing with snmc-seq2 hdf5 mapping, from fastq to allc, and generate QC stats.
plot/
    all plotting functions related to any step of analysis
tools/
    function for analysis and manipulating text based files (allc, bed, bigwig etc.).
local/
    unstable code that only used for local analysis, may be changed or removed in any time.
__main__.py
    enter point of yap

# TODO add api for useful functions, like scanpy
"""

__version__ = '0.1.2'
