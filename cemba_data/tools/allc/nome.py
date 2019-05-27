from .developing_map_to_region import map_to_sparse_chrom_bin
from .utilities import extract_allc_context
import subprocess
import glob


def _prepare_count_table(allc_path, out_prefix, chrom_size_file,
                         bin_size=500, mc_contexts=('GCYN', 'HCYN'),
                         remove_additional_chrom=False):
    extract_allc_context(allc_path=allc_path,
                         out_prefix=out_prefix,
                         mc_contexts=mc_contexts)

    map_to_sparse_chrom_bin(allc_path=allc_path,
                            out_prefix=out_prefix,
                            chrom_size_file=chrom_size_file,
                            remove_additional_chrom=remove_additional_chrom,
                            bin_size=bin_size)
    for path in glob.glob(out_prefix+'.extract*.tsv.gz*'):
        subprocess.run(['rm', 'f', path])
    return
