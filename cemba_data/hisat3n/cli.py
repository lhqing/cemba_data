import click

from .hisat3n_m3c import remove_overlap_read_parts


@click.command('remove_overlap_read_parts')
@click.argument('in_bam_path')
@click.argument('out_bam_path')
def _remove_overlap_read_parts(in_bam_path, out_bam_path):
    remove_overlap_read_parts(in_bam_path, out_bam_path)
    return


@click.group()
def _main():
    return


def main():
    _main.add_command(_remove_overlap_read_parts)
    _main()
    return
