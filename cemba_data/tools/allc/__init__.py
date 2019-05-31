"""
ALLC format related functions
"""
from .allc_to_bigwig import allc_to_bigwig
from .bam_to_allc import call_methylated_sites, batch_call_methylated_sites
from .merge_allc import merge_allc_files
from ._open import open_allc, open_bam, open_gz
from .utilities import standardize_allc, profile_allc, map_to_region, tabix_allc
from .extract_allc import extract_allc
