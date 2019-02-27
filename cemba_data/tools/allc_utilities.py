"""
This file is modified from methylpy https://github.com/yupenghe/methylpy

Author: Yupeng He



                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "{}"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright {2017} {Matthew D. Schultz and Yupeng He}

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import os
import shlex
import numpy as np
import subprocess
import multiprocessing
import gzip
import glob
import resource
import psutil
from .utilities import parse_mc_pattern, parse_chrom_size, genome_region_chunks, parse_file_paths
from concurrent.futures import ProcessPoolExecutor, as_completed
from .open import open_allc
import logging
import gc

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

# get the system soft and hard limit of file handle
SOFT, HARD = resource.getrlimit(resource.RLIMIT_NOFILE)
DEFAULT_MAX_ALLC = 150
PROCESS = psutil.Process(os.getpid())


def _increase_soft_fd_limit():
    """
    Increase soft file descriptor limit to hard limit,
    this is the maximum a process can do
    Use this in merge_allc, because for single cell, a lot of file need to be opened

    Some useful discussion
    https://unix.stackexchange.com/questions/36841/why-is-number-of-open-files-limited-in-linux
    https://docs.python.org/3.6/library/resource.html
    https://stackoverflow.com/questions/6774724/why-python-has-limit-for-count-of-file-handles/6776345
    """
    resource.setrlimit(resource.RLIMIT_NOFILE, (HARD, HARD))


def convert_allc_to_bigwig(input_allc_file,
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
    with _open_allc_file(input_allc_file) as f:
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


def merge_allc_files(allc_paths, out_path, chrom_size_file, bin_length=1000000, cpu=10):
    allc_files = parse_file_paths(allc_paths)

    try:
        with open(out_path, 'w'):
            pass
    except IOError:
        log.info("Can't create out_path")

    try:
        _check_tabix(allc_files)
        _batch_merge_allc_files_tabix(allc_files=allc_files,
                                      out_file=out_path,
                                      chrom_size_file=chrom_size_file,
                                      bin_length=bin_length,
                                      cpu=cpu)
    except FileNotFoundError:
        log.info('Some ALLC file do not have tabix, use the old merge function with .idx.'
                 'This function is slower and introduce much higher IO. Consider use bgzip/tabix')
        _merge_allc_files_idx(allc_files,
                              out_path,
                              num_procs=cpu,
                              mini_batch=DEFAULT_MAX_ALLC,
                              compress_output=True,
                              skip_snp_info=True)
    return


def _batch_merge_allc_files_tabix(allc_files, out_file, chrom_size_file, bin_length, cpu=10):
    regions = genome_region_chunks(chrom_size_file, bin_length=bin_length)
    log.info(f'Merge ALLC files with {cpu} processes')
    log.info(f'Split genome into {len(regions)} regions, each is {bin_length}bp')
    log.info(f'{len(allc_files)} to merge, the default ALLC file handel in 1 run is {DEFAULT_MAX_ALLC}')
    log.info(f'Process FH soft limit {SOFT}, hard limit {HARD}')

    _increase_soft_fd_limit()
    if len(allc_files) > DEFAULT_MAX_ALLC:
        # deal with too many allc files
        # merge in batches
        allc_fold = len(allc_files) // DEFAULT_MAX_ALLC + 1
        batch_allc = len(allc_files) // allc_fold + 1

        allc_files_batches = []
        out_paths = []
        for batch_id, i in enumerate(range(0, len(allc_files), batch_allc)):
            allc_files_batches.append(allc_files[i:i + batch_allc])
            out_paths.append(out_file + f'batch_{batch_id}.tmp.tsv.gz')
            if batch_id > 0:
                # merge last batch's merged allc into next batch
                allc_files_batches[batch_id].append(out_paths[batch_id - 1])
        out_paths[-1] = out_file
    else:
        allc_files_batches = [allc_files]
        out_paths = [out_file]

    log.info(f'Run merge in {len(allc_files_batches)} batches')
    log.info(' '.join(out_paths))

    for batch_num, (allc_files, out_file) in enumerate(zip(allc_files_batches, out_paths)):
        log.info(f'Run batch {batch_num}, '
                 f'{len(allc_files)} allc files, '
                 f'output to {out_file}')
        with open_allc(out_file, 'w', threads=3) as out_handle:

            # as_complete don't release, run total regions in sections to prevent too large memory
            parallel_section = 100
            for i in range(0, len(regions), parallel_section):
                cur_regions = regions[i:min(i + parallel_section, len(regions))]
                print(f'Running region from {cur_regions[0]} to {cur_regions[-1]}')
                with ProcessPoolExecutor(max_workers=cpu) as executor:
                    future_merge_result = {executor.submit(_merge_allc_files_tabix,
                                                           allc_files=allc_files,
                                                           out_file=None,
                                                           chrom_size_file=chrom_size_file,
                                                           query_region=region,
                                                           buffer_line_number=100000): region_id
                                           for region_id, region in enumerate(cur_regions)}
                    cur_id = 0
                    temp_dict = {}
                    # future may return in any order
                    # save future.result into temp_dict 1st
                    # write data in order by region_id
                    # so the out file is ordered
                    for future in as_completed(future_merge_result):
                        region_id = future_merge_result[future]
                        try:
                            temp_dict[region_id] = future.result()
                        except Exception as exc:
                            log.info('%r generated an exception: %s' % (region_id, exc))
                        else:
                            try:
                                while len(temp_dict) > 0:
                                    data = temp_dict.pop(cur_id)
                                    out_handle.write(data)
                                    log.info(f'write {cur_regions[cur_id]} Cached: {len(temp_dict)}, '
                                             f'Current memory size: {PROCESS.memory_info().rss/(1024**3):.2f}')
                                    cur_id += 1
                            except KeyError:
                                continue
                    # write last pieces of data
                    while len(temp_dict) > 0:
                        data = temp_dict.pop(cur_id)
                        out_handle.write(data)
                        log.info(f'write {regions[cur_id]} Cached: {len(temp_dict)}, '
                                 f'Current memory size: {PROCESS.memory_info().rss/(1024**3):.2f}')
                        cur_id += 1
                gc.collect()
        # after merge, tabix output
        log.info('Tabix output ALLC file')
        subprocess.run(['tabix', '-b', '2', '-e', '2', '-s', '1', out_file])
        log.info(f'Current memory size: {PROCESS.memory_info().rss/(1024**3):.2f}')
    log.info('Merge finished.')

    # remove all batch allc but the last (final allc)
    for out_file in out_paths[:-1]:  # last file is the final merged allc
        subprocess.run(shlex.split(f'rm -f {out_file} {out_file}.tbi'))
    return


def _check_tabix(allc_files):
    for allc_path in allc_files:
        if not os.path.exists(allc_path + '.tbi'):
            raise FileNotFoundError(f'Tabix for {allc_path} not found')
    return


def _merge_allc_files_tabix(allc_files,
                            out_file,
                            chrom_size_file,
                            query_region=None,
                            buffer_line_number=10000):
    # only use bgzip and tabix
    # automatically take care the too many file open issue
    # do merge iteratively if file number exceed limit
    # parallel in chrom_bin level, not chrom level

    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")
    chrom_size_dict = parse_chrom_size(chrom_size_file)
    all_chroms = list(chrom_size_dict.keys())
    if query_region is not None:
        if not isinstance(query_region, str):
            exit("query_region must be str or None")
        region_chroms = set([region.split(':')[0] for region in query_region.split(' ')])
        all_chroms = [chrom for chrom in all_chroms if chrom in region_chroms]
    processing_chrom = all_chroms.pop(0)

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    file_handles = [open_allc(allc_file,
                              region=query_region)
                    for allc_file in allc_files]
    if out_file is not None:
        out_handle = open_allc(out_file, 'w')
    else:
        out_handle = ''

    # merge allc files
    out = ""
    cur_chrom = ['NOT_A_CHROM' for _ in range(len(allc_files))]
    cur_pos = np.array([np.nan for _ in range(len(allc_files))])
    cur_fields = [None for _ in range(len(allc_files))]
    file_reading = np.array([True for _ in range(len(allc_files))])

    # init
    for index, allc_file in enumerate(allc_files):
        line = file_handles[index].readline()
        if line:
            fields = line.split("\t")
            cur_chrom[index] = fields[0]
            cur_pos[index] = int(fields[1])
            cur_fields[index] = fields
        else:
            # file handle read nothing, the file is empty
            file_reading[index] = False

    active_handle = np.array([True if chrom == processing_chrom else False
                              for chrom in cur_chrom])

    # merge
    line_count = 0
    while file_reading.sum() > 0:
        mc, cov = 0, 0
        genome_info = None
        # select index whose cur_pos is smallest among all active handle
        for index in np.where((cur_pos == np.nanmin(cur_pos[active_handle]))
                              & active_handle)[0]:
            mc += int(cur_fields[index][4])
            cov += int(cur_fields[index][5])
            if genome_info is None:
                genome_info = cur_fields[index][:4]

            # update
            line = file_handles[index].readline()
            if line:
                fields = line.split("\t")
                # judge if chrom changed between two lines
            else:
                # read to the end of a file
                fields = ['NOT_A_CHROM', 9999999999]
                file_reading[index] = False

            this_chrom = cur_fields[index][0]
            next_chrom = fields[0]
            if next_chrom == this_chrom:
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields
            else:
                # read to next chrom
                # handle became inactive
                active_handle[index] = False
                cur_chrom[index] = next_chrom
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields

                # if all handle became inactive, move processing_chrom to next
                if sum(active_handle) == 0:
                    if len(all_chroms) == 0:
                        break
                    processing_chrom = all_chroms.pop(0)
                    # and re-judge active handle
                    active_handle = np.array([True if chrom == processing_chrom else False
                                              for chrom in cur_chrom])

        # output
        out += '\t'.join(genome_info) + f'\t{mc}\t{cov}\t1\n'
        line_count += 1
        if line_count > buffer_line_number:
            if isinstance(out_handle, str):
                out_handle += out
            else:
                out_handle.write(out)
            line_count = 0
            out = ""
    # the last out
    for file_handle in file_handles:
        file_handle.close()
    if isinstance(out_handle, str):
        out_handle += out
        return out_handle
    else:
        out_handle.write(out)
        out_handle.close()
        return


def _merge_allc_files_idx(allc_files,
                          output_file,
                          num_procs=1,
                          mini_batch=DEFAULT_MAX_ALLC,
                          compress_output=True,
                          skip_snp_info=True):
    _increase_soft_fd_limit()
    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")
    # add .gz suffix
    if output_file[-3:] != ".gz":
        output_file += ".gz"

    # check index
    _index_allc_file_batch(allc_files,
                           cpu=min(num_procs, 10))
    log.info("Start merging")
    if not (num_procs > 1):
        _merge_allc_files_idx_minibatch(allc_files,
                                        output_file,
                                        query_chroms=None,
                                        mini_batch=mini_batch,
                                        compress_output=compress_output,
                                        skip_snp_info=skip_snp_info)
        _index_allc_file(output_file)
        return 0

    # parallel merging
    log.info("Getting chromosome names")
    chroms = set([])
    for allc_file in allc_files:
        c_p = _read_allc_index(allc_file)
        for chrom in c_p.keys():
            chroms.add(chrom)
    try:
        pool = multiprocessing.Pool(min(num_procs,
                                        len(chroms)))  # ,
        # int(float(mini_batch)/float(len(allc_files)))))
        log.info("Merging allc files")
        for chrom in chroms:
            pool.apply_async(_merge_allc_files_idx_minibatch,
                             (),
                             {"allc_files": allc_files,
                              "output_file": output_file + "_" + str(chrom) + ".tsv",
                              "query_chroms": chrom,
                              "mini_batch": mini_batch,
                              "compress_output": False,
                              "skip_snp_info": skip_snp_info,
                              })
        pool.close()
        pool.join()
        # output
        log.info("Merging outputs")
        file_list = []
        for chrom in chroms:
            file_list.append(output_file + "_" + chrom + ".tsv")
        cat_process = subprocess.Popen(['cat'] + file_list,
                                       stdout=subprocess.PIPE)
        with open(output_file, 'w') as f:
            subprocess.run(['pigz'], stdin=cat_process.stdout,
                           stdout=f)
    except:
        log.info("Failed to merge using multiple processors. " +
                 "Do minibatch merging using single processor.")
        _merge_allc_files_idx_minibatch(allc_files,
                                        output_file,
                                        mini_batch=mini_batch,
                                        compress_output=compress_output,
                                        skip_snp_info=skip_snp_info)
    # remove temporary files
    for chrom in set(chroms):
        tmp_file = glob.glob(output_file + "_" + str(chrom) + ".tsv")
        if tmp_file:
            tmp_file = tmp_file[0]
            subprocess.run(["rm", tmp_file])
    # index output allc file
    _index_allc_file(output_file)
    return 0


def _merge_allc_files_idx_minibatch(allc_files,
                                    output_file,
                                    query_chroms=None,
                                    mini_batch=100,
                                    compress_output=False,
                                    skip_snp_info=True):
    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")
    # merge all files at once
    try:
        _merge_allc_files_idx_worker(allc_files=allc_files,
                                     output_file=output_file,
                                     query_chroms=query_chroms,
                                     compress_output=compress_output,
                                     skip_snp_info=skip_snp_info)
        return 0
    except:
        log.info("Failed to merge all allc files at once. Do minibatch merging")

    # init
    remaining_allc_files = list(allc_files[mini_batch:])
    output_tmp_file = output_file + ".tmp"
    _merge_allc_files_idx_worker(allc_files=allc_files[:mini_batch],
                                 output_file=output_file,
                                 query_chroms=query_chroms,
                                 compress_output=compress_output,
                                 skip_snp_info=skip_snp_info)
    # batch merge
    while len(remaining_allc_files) > 0:
        processing_allc_files = [output_file]
        while len(remaining_allc_files) > 0 \
                and len(processing_allc_files) < mini_batch:
            processing_allc_files.append(remaining_allc_files.pop())
        _merge_allc_files_idx_worker(allc_files=processing_allc_files,
                                     output_file=output_tmp_file,
                                     query_chroms=query_chroms,
                                     compress_output=compress_output,
                                     skip_snp_info=skip_snp_info)
        subprocess.check_call(["mv", output_tmp_file, output_file])
        _index_allc_file(output_file)
    return 0


def _merge_allc_files_idx_worker(allc_files,
                                 output_file,
                                 query_chroms=None,
                                 compress_output=False,
                                 skip_snp_info=True,
                                 buffer_line_number=100000):
    # User input checks
    if not isinstance(allc_files, list):
        exit("allc_files must be a list of string(s)")

    # scan allc file to set up a table for fast look-up of lines belong
    # to different chromosomes
    fhandles = []
    chrom_pointer = []
    chroms = set([])
    try:
        for index, allc_file in enumerate(allc_files):
            fhandles.append(_open_allc_file(allc_file))
            chrom_pointer.append(_read_allc_index(allc_file))
            for chrom in chrom_pointer[index].keys():
                chroms.add(chrom)
    except:
        for f in fhandles:
            f.close()
        raise
        # exit() # exit due to failure of openning all allc files at once
    if query_chroms is not None:
        if isinstance(query_chroms, list):
            chroms = query_chroms
        else:
            chroms = [query_chroms]
    chroms = list(map(str, chroms))
    # output
    if compress_output:
        g = gzip.open(output_file, 'wt')
    else:
        g = open(output_file, 'w')

    # merge allc files
    out = ""
    for chrom in chroms:
        cur_pos = np.array([np.nan for index in range(len(allc_files))])
        cur_fields = []
        num_remaining_allc = 0
        # init
        for index, allc_file in enumerate(allc_files):
            fhandles[index].seek(chrom_pointer[index].get(chrom, 0))
            line = fhandles[index].readline()
            fields = line.split("\t")
            cur_fields.append(None)
            if fields[0] == chrom:
                cur_pos[index] = int(fields[1])
                cur_fields[index] = fields
                num_remaining_allc += 1

        # check consistency of SNP information
        context_len = None
        if not skip_snp_info:
            for index, allc_file in enumerate(allc_files):
                # skip allc with no data in this chromosome
                if cur_fields[index] is None:
                    continue
                # SNP information is missing
                if len(cur_fields[index]) < 9:
                    skip_snp_info = True
                    break
                # init contex_len
                if context_len is None:
                    context_len = len(cur_fields[index][7].split(","))
                # check whether contex length are the same
                if context_len != len(cur_fields[index][7].split(",")) or \
                        context_len != len(cur_fields[index][8].split(",")):
                    log.info("Inconsistent sequence context length: "
                             + allc_file + "\n")
        # merge
        line_counts = 0
        while num_remaining_allc > 0:
            mc, h = 0, 0
            if not skip_snp_info:
                matches = [0 for _ in range(context_len)]
                mismatches = list(matches)
            c_info = None
            for index in np.where(cur_pos == np.nanmin(cur_pos))[0]:
                mc += int(cur_fields[index][4])
                h += int(cur_fields[index][5])
                if not skip_snp_info:
                    for ind, match, mismatch in zip(range(context_len),
                                                    cur_fields[index][7].split(","),
                                                    cur_fields[index][8].split(",")):
                        matches[ind] += int(match)
                        mismatches[ind] += int(mismatch)
                if c_info is None:
                    c_info = "\t".join(cur_fields[index][:4])
                # update
                line = fhandles[index].readline()
                fields = line.split("\t")
                if fields[0] == chrom:
                    cur_pos[index] = int(fields[1])
                    cur_fields[index] = fields
                else:
                    cur_pos[index] = np.nan
                    num_remaining_allc -= 1
            # output
            if not skip_snp_info:
                out += c_info + "\t" + str(mc) + "\t" + str(h) + "\t1" + "\t" \
                       + ",".join(map(str, matches)) + "\t" + ",".join(map(str, mismatches)) + "\n"
            else:
                out += c_info + "\t" + str(mc) + "\t" + str(h) + "\t1\n"
            line_counts += 1
            if line_counts > buffer_line_number:
                g.write(out)
                line_counts = 0
                out = ""
    if line_counts > 0:
        g.write(out)
    g.close()
    for index in range(len(allc_files)):
        fhandles[index].close()
    return 0


def _index_allc_file_batch(allc_files, cpu=1):
    if isinstance(allc_files, str):
        allc_files = glob.glob(allc_files)

    if cpu == 1:
        for allc_file in allc_files:
            _index_allc_file(allc_file)
    else:
        pool = multiprocessing.Pool(min(cpu, len(allc_files), 48))
        for allc_file in allc_files:
            pool.apply_async(_index_allc_file, (allc_file,))
        pool.close()
        pool.join()
    return 0


def _index_allc_file(allc_file):
    index_file = str(allc_file) + '.idx'
    if os.path.exists(index_file):
        # backward compatibility
        with open(index_file) as f:
            last_line = f.readlines()[-1]
            if last_line == "#eof\n":
                return

    if allc_file.endswith('gz'):  # works for .gz, .bgz
        f = subprocess.Popen(['zcat', allc_file],
                             stdout=subprocess.PIPE,
                             encoding='utf8').stdout
    else:
        f = open(allc_file)

    index_lines = []
    cur_chrom = "TOTALLY_NOT_A_CHROM"
    cur_start = cur_chrom + '\t'
    cur_pointer = 0
    # check header
    first_line = True
    for line in f:
        if first_line:
            fields = line.split("\t")
            first_line = False
            try:
                int(fields[1])
                int(fields[4])
                int(fields[5])
                # no header, continue to start from the beginning of allc file
            except ValueError:
                # find header, skip it
                cur_pointer += len(line)
                continue
        if not line.startswith(cur_start):
            fields = line.split("\t")
            index_lines.append(fields[0] + "\t" + str(cur_pointer) + "\n")
            cur_chrom = fields[0]
            cur_start = cur_chrom + '\t'
        cur_pointer += len(line)
    # backward compatibility
    index_lines.append("#eof\n")
    f.close()
    with open(index_file, 'w') as idx:
        idx.writelines(index_lines)
    return


def _read_allc_index(allc_file):
    index_file = str(allc_file) + ".idx"
    _index_allc_file(allc_file)
    with open(index_file, 'r') as f:
        chrom_pointer = {}
        for line in f:
            if line[0] != '#':
                fields = line.rstrip().split("\t")
                chrom_pointer[fields[0]] = int(fields[1])
    return chrom_pointer


def _open_allc_file(allc_file):
    if allc_file.endswith(".gz"):
        f = gzip.open(allc_file, 'rt')
    else:
        f = open(allc_file, 'r')
    return f
