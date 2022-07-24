import pysam
import pathlib
import cemba_data
import subprocess
from ..utilities import get_configuration


def bam_read_to_fastq_read(read):
    if read.is_read1:
        read_type = '1'
    else:
        read_type = '2'

    fastq_record = f"@{read.qname}_{read_type}\n" \
                   f"{read.query_sequence}\n" \
                   f"+\n" \
                   f"{read.qual}\n"
    return fastq_record


def separate_unique_and_multi_align_reads(in_bam_path,
                                          out_unique_path,
                                          out_multi_path,
                                          out_unmappable_path=None,
                                          unmappable_format='auto',
                                          mapq_cutoff=10,
                                          qlen_cutoff=30):
    """
    Separate unique aligned, multi-aligned, and unaligned reads from hisat-3n bam file.

    Parameters
    ----------
    in_bam_path
        Path to hisat-3n bam file.
    out_unique_path
        Path to output unique aligned bam file.
    out_multi_path
        Path to output multi-aligned bam file.
    out_unmappable_path
        Path to output unmappable file.
    unmappable_format
        Format of unmappable file, only "bam" and "fastq" supported.
    mapq_cutoff
        MAPQ cutoff for uniquely aligned reads,
        note that for hisat-3n, unique aligned reads always have MAPQ=60
    qlen_cutoff
        read length cutoff for any reads
    Returns
    -------
    None
    """
    if out_unmappable_path is not None:
        if unmappable_format == 'auto':
            if out_unmappable_path.endswith('.bam'):
                unmappable_format = 'bam'
            elif out_unmappable_path.endswith('.fastq'):
                unmappable_format = 'fastq'
            else:
                raise ValueError(f'Unmappable format {unmappable_format} not supported.')
        else:
            if unmappable_format not in ['bam', 'fastq']:
                raise ValueError(f'Unmappable format {unmappable_format} not supported.')

    with pysam.AlignmentFile(in_bam_path, index_filename=None) as bam:
        header = bam.header
        with pysam.AlignmentFile(out_unique_path, header=header, mode='wb') as unique_bam, \
                pysam.AlignmentFile(out_multi_path, header=header, mode='wb') as multi_bam:
            if out_unmappable_path is not None:
                if unmappable_format == 'bam':
                    unmappable_file = pysam.AlignmentFile(out_unmappable_path, header=header, mode='wb')
                else:
                    unmappable_file = open(out_unmappable_path, 'w')
            else:
                unmappable_file = None

            for read in bam:
                # skip reads that are too short
                if read.qlen < qlen_cutoff:
                    continue

                if read.mapq > mapq_cutoff:
                    unique_bam.write(read)
                elif read.mapq > 0:
                    multi_bam.write(read)
                else:
                    # unmappable reads
                    if unmappable_file is not None:
                        if unmappable_format == 'bam':
                            unmappable_file.write(read)
                        else:
                            unmappable_file.write(bam_read_to_fastq_read(read))

            if unmappable_file is not None:
                unmappable_file.close()
    return


def convert_hisat_bam_strandness(in_bam_path, out_bam_path):
    with pysam.AlignmentFile(in_bam_path) as in_bam, \
            pysam.AlignmentFile(out_bam_path, header=in_bam.header, mode='wb') as out_bam:
        for read in in_bam:
            if read.get_tag('YZ') == '+':
                read.is_forward = True
                if read.is_paired:
                    read.mate_is_forward = True
            else:
                read.is_forward = False
                if read.is_paired:
                    read.mate_is_forward = False
            out_bam.write(read)
    return


def make_snakefile_hisat3n(output_dir):
    output_dir = pathlib.Path(output_dir)

    mapping_config_name = list(output_dir.glob('mapping_config.*'))[0].name

    config = get_configuration(output_dir / mapping_config_name)
    try:
        mode = config['mode']
    except KeyError:
        raise KeyError('mode not found in the config file.')

    skip_dirs = ['stats', 'snakemake', 'scool']
    mapping_job_dirs = [p for p in output_dir.glob('*')
                        if p.is_dir() and (p.name not in skip_dirs)]

    snakemake_dir = output_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True)
    stats_dir = output_dir / 'stats'
    stats_dir.mkdir(exist_ok=True)

    package_dir = cemba_data.__path__[0]
    snakefile_path = f'{package_dir}/hisat3n/snakefile/{mode.lower()}.smk'
    if not pathlib.Path(snakefile_path).exists():
        print('Possible snakefile templates:')
        for p in pathlib.Path(f'{package_dir}/hisat3n/snakefile/').glob('*.smk'):
            print(p)
        raise ValueError(f'Mode {mode} not supported, '
                         f'because Snakefile {snakefile_path} not found.')

    for p in mapping_job_dirs:
        subprocess.run(['cp', f'{output_dir}/{mapping_config_name}',
                        f'{p}/{mapping_config_name}'], check=True)
        subprocess.run(['cp', snakefile_path, f'{p}/Snakefile'], check=True)

    # leave a flag to indicate using hisat-3n pipeline
    subprocess.run(['touch', f'{output_dir}/snakemake/hisat3n'], check=True)
    return
