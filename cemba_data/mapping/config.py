import pathlib

import cemba_data
from ..utilities import MAPPING_MODE_CHOICES

# Load defaults
PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])


def print_default_mapping_config(mode,
                                 barcode_version,
                                 genome_fasta,
                                 bismark_ref=None,
                                 hisat3n_dna_ref=None,
                                 hisat3n_rna_ref=None,
                                 star_ref=None,
                                 gtf=None,
                                 nome=False,
                                 chrom_size_path=None):
    mode = mode.lower()
    if mode not in MAPPING_MODE_CHOICES:
        raise ValueError(f'Unknown mode {mode}')

    barcode_version = barcode_version.upper()
    if barcode_version not in ['V1', 'V2']:
        raise ValueError(f'Unknown mode {barcode_version}')

    if bismark_ref is not None:
        bismark_ref = pathlib.Path(bismark_ref).absolute()
    else:
        if hisat3n_rna_ref is None or hisat3n_dna_ref is None:
            raise ValueError('bismark_ref is required if hisat3n_rna_ref and hisat3n_dna_ref are not specified.')
        hisat3n_rna_ref = pathlib.Path(hisat3n_rna_ref).absolute()
        hisat3n_dna_ref = pathlib.Path(hisat3n_dna_ref).absolute()

    if mode == 'mct':
        if star_ref is None:
            if hisat3n_rna_ref is None:
                raise ValueError('star_ref or hisat3n_rna_ref is required if mode is mct.')
        else:
            star_ref = pathlib.Path(star_ref).absolute()
        if gtf is None:
            raise ValueError('gtf must be provided when mode is mct.')
        gtf = pathlib.Path(gtf).absolute()

    if chrom_size_path is None:
        raise ValueError('chrom_size_path must be provided.')
    chrom_size_path = pathlib.Path(chrom_size_path).absolute()

    if mode == 'm3c':
        pass

    if mode == '4m':
        if (star_ref is None) and (hisat3n_rna_ref is None):
            raise ValueError('star_ref or hisat3n_rna_ref is required if mode is mct.')
        star_ref = pathlib.Path(star_ref).absolute()

        if gtf is None:
            raise ValueError('gtf must be provided when mode is mct.')
        gtf = pathlib.Path(gtf).absolute()

    genome_fasta = pathlib.Path(genome_fasta).absolute()

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
        if hisat3n_rna_ref is None:
            config_content = config_content.replace('CHANGE_THIS_TO_YOUR_STAR_REFERENCE_DIR', str(star_ref))
        else:
            config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_RNA_REFERENCE',
                                                    str(hisat3n_rna_ref))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_GENE_ANNOTATION_GTF', str(gtf))
    elif mode == 'm3c':
        config_path = PACKAGE_DIR / 'files/default_config/mapping_config_m3c.ini'
        with open(config_path) as f:
            config_content = f.read()
    elif mode == '4m':
        config_path = PACKAGE_DIR / 'files/default_config/mapping_config_4m.ini'
        with open(config_path) as f:
            config_content = f.read()
        if hisat3n_rna_ref is None:
            config_content = config_content.replace('CHANGE_THIS_TO_YOUR_STAR_REFERENCE_DIR', str(star_ref))
        else:
            config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_RNA_REFERENCE',
                                                    str(hisat3n_rna_ref))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_GENE_ANNOTATION_GTF', str(gtf))
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_CHROM_SIZE_PATH', str(chrom_size_path))
    else:
        raise

    config_content = config_content.replace('CHANGE_THIS_TO_YOUR_CHROM_SIZE_PATH', str(chrom_size_path))
    config_content = config_content.replace('USE_CORRECT_BARCODE_VERSION_HERE', barcode_version)
    if hisat3n_dna_ref is None:
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_BISMARK_REFERENCE_DIR', str(bismark_ref))
    else:
        config_content = config_content.replace('CHANGE_THIS_TO_YOUR_HISAT3N_DNA_REFERENCE', str(hisat3n_dna_ref))
    config_content = config_content.replace('CHANGE_THIS_TO_YOUR_REFERENCE_FASTA', str(genome_fasta))
    print(config_content)
    return
