import pysam
import pathlib
import cemba_data
import subprocess
from ..utilities import get_configuration


def separate_unique_and_multi_align_reads(in_bam_path,
                                          out_unique_path,
                                          out_multi_path,
                                          out_unmappable_path=None,
                                          mapq_cutoff=10,
                                          qlen_cutoff=30):
    # TODO: make sure how to separate multi-align reads? or reads map to repeat regions in the genome?
    # TODO right now, we are just using mapq == 1 as multi-align reads, but this might not be right

    with pysam.AlignmentFile(in_bam_path, index_filename=None) as bam:
        header = bam.header
        with pysam.AlignmentFile(out_unique_path, header=header, mode='wb') as unique_bam, \
                pysam.AlignmentFile(out_multi_path, header=header, mode='wb') as multi_bam:
            if out_unmappable_path is not None:
                unmappable_bam = pysam.AlignmentFile(out_unmappable_path, header=header, mode='wb')
            else:
                unmappable_bam = None

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
                    if unmappable_bam is not None:
                        unmappable_bam.write(read)

            if unmappable_bam is not None:
                unmappable_bam.close()
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
    mapping_job_dirs = [p for p in output_dir.glob('*') if p.is_dir() and (p.name not in skip_dirs)]

    snakemake_dir = output_dir / 'snakemake'
    snakemake_dir.mkdir(exist_ok=True)
    stats_dir = output_dir / 'stats'
    stats_dir.mkdir(exist_ok=True)

    package_dir = cemba_data.__path__[0]
    snakefile_path = f'{package_dir}/hisat3n/snakefile/{mode.lower()}.Snakefile'
    if not pathlib.Path(snakefile_path).exists():
        print('Possible snakefile templates:')
        for p in pathlib.Path(f'{package_dir}/hisat3n/snakefile/').glob('Snakefile.*'):
            print(p)
        raise ValueError(f'Mode {mode} not supported, because Snakefile {snakefile_path} not found.')

    for p in mapping_job_dirs:
        subprocess.run(['cp', f'{output_dir}/{mapping_config_name}', f'{p}/{mapping_config_name}'], check=True)
        subprocess.run(['cp', snakefile_path, f'{p}/Snakefile'], check=True)

    # leave a flag to indicate using hisat-3n pipeline
    subprocess.run(['touch', f'{output_dir}/snakemake/hisat3n'], check=True)
    return
