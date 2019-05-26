"""
This file is modified from methylpy https://github.com/yupenghe/methylpy

Original author: Yupeng He
"""

import shlex
import subprocess
from ..utilities import parse_mc_pattern
from .open import open_allc
import logging

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())


def _allc_to_bigwig(input_allc_file,
                    output_file,
                    chrom_size_file,
                    mc_type="CGN",
                    bin_size=100,
                    path_to_wigtobigwig="",
                    min_bin_sites=0,
                    min_bin_cov=0,
                    max_site_cov=None,
                    min_site_cov=0,
                    add_chr_prefix=False):
    # TODO add option for only convert cov bigwig
    if not isinstance(mc_type, list):
        if isinstance(mc_type, str):
            mc_type = [mc_type]
        else:
            exit("mc_type must be a list of string(s)")
    # TODO add supports for multiple mc_type
    mc_type = mc_type[0]

    if len(path_to_wigtobigwig):
        path_to_wigtobigwig += "/"

    mc_class = parse_mc_pattern(mc_type)

    chrom_end = {}

    # chromosome size
    f = open(chrom_size_file, 'r')
    g = open(output_file + ".chrom_size", 'w')
    for line in f:
        fields = line.split("\t")
        if not fields[0].startswith("chr"):
            fields[0] = "chr" + fields[0]
        chrom_end[fields[0]] = int(fields[1])
        g.write(fields[0] + "\t" + fields[1] + "\n")
    g.close()

    # prepare wig file
    cur_chrom = ""
    bin_start, bin_end = 0, 0
    bin_mc, bin_h, bin_site = 0, 0, 0
    g = open(output_file + ".wig", 'w')
    with open_allc(input_allc_file) as f:
        for line in f:
            fields = line.split("\t")
            if fields[3] not in mc_class:
                continue
            pos = int(fields[1]) - 1
            if cur_chrom != fields[0] or pos >= bin_end:
                try:
                    if fields[0].startswith("chr"):
                        cur_chrom_end = chrom_end[fields[0]]
                    else:
                        cur_chrom_end = chrom_end["chr" + fields[0]]
                except KeyError:
                    # chrom not in chrom size file
                    continue

                if bin_h > 0 and bin_site >= min_bin_sites and bin_h >= min_bin_cov:
                    mc_level = str(float(bin_mc) / float(bin_h))
                    if add_chr_prefix and not cur_chrom.startswith("chr"):
                        g.write("\t".join(["chr" + cur_chrom,
                                           str(bin_start),
                                           str(bin_end),
                                           mc_level]) + "\n")
                    else:
                        g.write("\t".join([cur_chrom,
                                           str(bin_start),
                                           str(bin_end),
                                           mc_level]) + "\n")

                if pos >= cur_chrom_end:
                    cur_chrom = fields[0]
                    bin_mc, bin_h, bin_site = 0, 0, 0
                    continue
                # reset
                cur_chrom = fields[0]
                bin_mc, bin_h, bin_site = 0, 0, 0
                bin_end = int(float(pos) / float(bin_size) + 1) * bin_size
                bin_start = bin_end - bin_size
                if bin_end > cur_chrom_end:
                    bin_end = cur_chrom_end
            # update mc, h and site
            h_site = int(fields[5])
            if h_site >= min_site_cov \
                    and (max_site_cov is None or h_site <= max_site_cov):
                bin_mc += int(fields[4])
                bin_h += h_site
                bin_site += 1
    if bin_h > 0 and bin_site >= min_bin_sites and bin_h >= min_bin_cov:
        mc_level = str(float(bin_mc) / float(bin_h))
        if add_chr_prefix and not cur_chrom.startswith("chr"):
            g.write("\t".join(["chr" + cur_chrom,
                               str(bin_start),
                               str(bin_end),
                               mc_level]) + "\n")
        else:
            g.write("\t".join([cur_chrom,
                               str(bin_start),
                               str(bin_end),
                               mc_level]) + "\n")
    g.close()

    # generate bigwig file
    subprocess.check_call(shlex.split(path_to_wigtobigwig + "wigToBigWig "
                                      + "%s.wig " % output_file
                                      + "%s.chrom_size " % output_file
                                      + output_file), stderr=subprocess.PIPE)
    subprocess.check_call(shlex.split("rm " + output_file + ".wig " + output_file + ".chrom_size"))


def allc_to_bigwig(allc_path, out_path, chrom_size, mc_type='CGN'):
    # TODO add allc to bigwig COV version, not calculate mC but only compute cov
    _allc_to_bigwig(allc_path,
                    out_path,
                    chrom_size,
                    mc_type=mc_type,
                    bin_size=100,
                    path_to_wigtobigwig="",
                    min_bin_sites=0,
                    min_bin_cov=0,
                    max_site_cov=None,
                    min_site_cov=0,
                    add_chr_prefix=True)
    return
