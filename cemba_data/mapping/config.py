import pathlib

import cemba_data
from ..utilities import MAPPING_MODE_CHOICES

# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def print_default_mapping_config(mode, barcode_version, bismark_ref, genome_fasta,
                                 star_ref=None, gtf=None, nome=False):
    mode = mode.lower()
    if mode not in MAPPING_MODE_CHOICES:
        raise ValueError(f'Unknown mode {mode}')

    barcode_version = barcode_version.upper()
    if barcode_version not in ['V1', 'V2']:
        raise ValueError(f'Unknown mode {barcode_version}')

    bismark_ref = pathlib.Path(bismark_ref).absolute()
    if not (bismark_ref / 'Bisulfite_Genome').exists():
        raise FileNotFoundError(f"{bismark_ref / 'Bisulfite_Genome'} not exist")

    if mode == 'mct':
        if (star_ref is None) or (gtf is None):
            raise ValueError('star_ref and gtf must be provided when mode is mct.')

        star_ref = pathlib.Path(star_ref).absolute()
        if not star_ref.exists():
            raise FileNotFoundError(f"{star_ref} not exist")

        gtf = pathlib.Path(gtf).absolute()
        if not gtf.exists():
            raise FileNotFoundError(f'{gtf} not exist')

    genome_fasta = pathlib.Path(genome_fasta).absolute()
    if not genome_fasta.exists():
        raise FileNotFoundError(f'{genome_fasta} not exist')

    if mode == 'mc':
        if nome:
            config_path = PACKAGE_DIR / 'files/default_config/mapping_config_nome.ini'
        else:
            config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mc.ini'
        with open(config_path) as f:
            config_content = f.read()
    elif mode == 'mct':
        if nome:
            config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mct-nome.ini'
        else:
            config_path = PACKAGE_DIR / 'files/default_config/mapping_config_mct.ini'
        with open(config_path) as f:
            config_content = f.read()
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_STAR_REFERENCE_DIR', str(star_ref))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_GENE_ANNOTATION_GTF', str(gtf))
    else:
        raise

    config_content = config_content.replace('USE_CORRECT_BARCODE_VERSION_HERE', barcode_version)
    config_content = config_content.replace('CHANGE_THIS_TO_YOUR_BISMARK_REFERENCE_DIR', str(bismark_ref))
    config_content = config_content.replace('CHANGE_THIS_TO_YOUR_REFERENCE_FASTA', str(genome_fasta))
    print(config_content)
    return
