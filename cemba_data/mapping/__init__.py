"""
Code order

1. plateinfo_and_samplesheet.py
Prepare plateinfo file and samplesheet file, this is where the name pattern start

2. fastq_dataframe.py
Prepare fastq_dataframe based on consistent pattern from samplesheet

3. demultiplex.py
Random index demultiplex to get single cell fastq files

4. merge_lane.py
QC raw fastq to get trimmed and filtered fastq

"""