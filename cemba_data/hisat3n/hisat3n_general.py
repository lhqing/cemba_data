import pysam


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
