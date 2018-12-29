"""
Test functions in cemba_data.mapping
Mapping related functions from FASTQ to ALLC, together with the mapping statistics

Testing related files should exist in ./_data/mapping
"""

from cemba_data.mapping.fastq import demultiplex, fastq_qc
from cemba_data.mapping.bam import bam_qc
from cemba_data.mapping.bismark import bismark
from cemba_data.mapping.allc import call_methylated_sites
from cemba_data.mapping.pipeline import pipeline, validate_fastq_dataframe, summary_pipeline_stat


def test_demultiplex():
    return


def test_fastq_qc():
    return


def test_bam_qc():
    return


def test_bismark():
    return


def test_call_methylated_sites():
    return


def test_pipeline():
    return


def test_validate_fastq_dataframe():
    # use the provided fastq_dataframe, make sure the dataframe validation function worked
    return


def test_summary_pipeline_stat():
    return
