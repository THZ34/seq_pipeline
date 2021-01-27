import multiprocessing
import math


def fastp(sample, kwargs):
    """Shell command template
    Quality control with fastp

    :param sample: sample name
    :param suffix: suffix
    :param trim: trim or not (default trim)
    :return:
    """
    # fastp quality control
    # contain adapter trim
    # output filtered fastq and report .json/.html
    suffix = kwargs['suffix']
    trim = kwargs['trim']
    if trim == False:
        trim = '-A'
    elif trim == True:
        trim = ''
    command = "fastp " \
              "-i fastqs/%s%s " \
              "-o fastqs/%s_qc%s " \
              "-I fastqs/%s%s " \
              "-O fastqs/%s_qc%s " \
              "-j fastqs/%s.json " \
              "-h fastqs/%s.html " \
              "--thread %s " \
              "%s" % (
                  sample, suffix[0], sample, suffix[0], sample, suffix[1], sample, suffix[1], sample, sample,
                  max(1, math.floor(multiprocessing.cpu_count() * 0.9)), trim)

    return command
