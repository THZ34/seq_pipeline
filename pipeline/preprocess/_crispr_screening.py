import pandas as pd
from typing import List, Tuple, TypedDict, AnyStr
import logging


def mageck_count(samples, suffix, library, sgrna_len, out, trim=None, control_sgrna=None,
                 sample_labels=None):
    """
    :param samples:
    :param suffix:
    :param library:
    :param sgrna_len:
    :param trim:
    :param out:
    :param control_sgrna:
    :return:
    """
    #
    filelist_1 = []
    filelist_2 = []
    sample_list = []
    for sample in samples:
        if type(sample) == str:
            filelist_1.append(sample + suffix[0])
            filelist_2.append(sample + suffix[1])
            sample_list.append(sample)
        elif type(sample) == tuple:
            sample_list.append(tuple[0])
            rep_filelist_1 = []
            rep_filelist_2 = []
            for rep in sample:
                rep_filelist_1.append(rep + suffix[0])
                rep_filelist_2.append(rep + suffix[1])
            filelist_1.append(','.join(rep_filelist_1))
            filelist_2.append(','.join(rep_filelist_2))

    if not sample_labels:
        sample_labels = sample_list

    fastq_1 = '--fastq ' + ' '.join(filelist_1)
    fastq_2 = '--fastq-2 ' + ' '.join(filelist_2)
    sample_labels = '--sample-label ' + ' '.join(sample_labels)
    library = '-l ' + library
    sgrna_len = '--sgrna-len ' + str(sgrna_len)
    out = '-n ' + out
    if trim:
        trim = '--trim ' + str(trim)
    else:
        trim = ''

    if control_sgrna:
        control_sgrna = '--control-sgrna ' + control_sgrna
    else:
        control_sgrna = ''
    command = 'mageck count %s %s %s %s %s %s %s ' % (
        fastq_1, fastq_2, sample_labels, library, sgrna_len, control_sgrna, out)

    logging.debug(filelist_1)
    logging.debug(filelist_2)
    logging.debug(sample_list)

    return command


def mageck_count_report():
    pass
