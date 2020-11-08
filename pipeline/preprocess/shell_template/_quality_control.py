def fastp(sample, kwargs):
    """Shell command template
    Quality control with fastp

    :param sample:
    :param suffix:
    :return:
    """
    # fastp quality control
    # contain adapter trim
    # output filtered fastq and report .json/.html
    suffix = kwargs['suffix']
    command = "fastp " \
              "-i fastq/%s%s " \
              "-I fastq/%s_qc%s " \
              "-o fastq/%s%s " \
              "-O fastq/%s_qc%s " \
              "-j fastp/%s.json " \
              "-h fastp/%s.html" % (
                  sample, suffix[0], sample, suffix[0], sample, suffix[1], sample, suffix[1], sample, sample)

    return command
