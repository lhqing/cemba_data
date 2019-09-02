import pathlib


def pipeline_mc_info(output_dir, config_path, mct=False, qsub=True):
    # create directory
    output_dir = pathlib.Path(output_dir)
    fastq_dir = output_dir / 'fastq'
    bam_dir = output_dir / 'star_bam'
    bam_dir.mkdir(exist_ok=True, parents=True)

    # prepare STAR mapping

    # runner STAR mapping

    # summary STAR mapping

    # prepare STAR BAM QC

    # runner STAR BAM QC

    # summary STAR BAM QC

    # Prepare select RNA reads

    # runner select RNA reads

    # summary select RNA reads

    # feature count

    # final summary
