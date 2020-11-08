from ._chip import *
from ._align import *
from ._quality_control import *
from ._crispr_screening import *
from .._shell_command_generater import shell
from ._rna import *


def crispr_screening(samples, suffix, library, sgrna_len=None, trim=None, control_sgrna=None,
                     sample_labels=None, design_matrix=None):
    """Shell command template

    :param samples:
    :param suffix:
    :param library:
    :param sgrna_len:
    :param out:
    :param trim:
    :param control_sgrna:
    :param sample_labels:
    :param design_matrix:
    :return:
    """
    command_list = []
    command_list.extend(shell(samples, fastp, suffix=suffix))
    command_list.extend(shell(samples, mageck_count, suffix=suffix, library=library, sgrna_len=sgrna_len, trim=trim,
                              control_sgrna=control_sgrna, sample_labels=sample_labels, design_matrix=design_matrix))
    return command_list


def cuttag(samples, suffix, reference):
    """Shell command template

    :param samples:
    :param suffix:
    :param reference:
    :return:
    """
    command_list = []
    command_list.extend(shell(samples, fastp, suffix=suffix))
    command_list.extend(shell(samples, bowtie2_cuttag, suffix=suffix, reference=reference))
    command_list.extend(shell(samples, macs2))
    return command_list


def rna_salmon(samples, suffix, reference):
    """"""
    command_list = []
    command_list.extend(shell(samples, fastp, suffix=suffix))
    command_list.extend(shell(samples, salmon,reference=reference))
