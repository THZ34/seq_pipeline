import logging


def mageck_count(samples, kwargs):
    """ Template of mageck count
    :param samples:
    :param suffix:
    :param library:
    :param sgrna_len: Length of sgRNA
    :param trim: Length for trim-5
    :param control_sgrna: Gene name of control sgRNA
    :return: mageck command
    """

    # args parse
    try:
        suffix = kwargs['suffix']
        library = kwargs['library']
    except:
        lost_arg = list({'suffix', 'library'} - set(kwargs.keys()))
        logging.error('Template required parameter missing :' + ' '.join(lost_arg))
        return 'Template required parameter missing'
    sgrna_len = kwargs.get('sgrna_len')
    trim = kwargs.get('trim')
    control_sgrna = kwargs.get('control_sgrna')
    sample_labels = kwargs.get('sample_labels')

    # process use quality controled fastq
    filelist_1 = []
    filelist_2 = []
    sample_list = []
    for sample in samples:
        if type(sample) == str:
            filelist_1.append(sample + '_qc' + suffix[0])
            filelist_2.append(sample + '_qc' + suffix[1])
            sample_list.append(sample)
        elif type(sample) == tuple:
            sample_list.append(sample[0])
            rep_filelist_1 = []
            rep_filelist_2 = []
            for rep in sample:
                rep_filelist_1.append(rep + '_qc' + suffix[0])
                rep_filelist_2.append(rep + '_qc' + suffix[1])
            filelist_1.append(','.join(rep_filelist_1))
            filelist_2.append(','.join(rep_filelist_2))

    if not sample_labels:
        sample_labels = sample_list

    fastq_1 = '--fastq ' + ' '.join(filelist_1)
    fastq_2 = '--fastq-2 ' + ' '.join(filelist_2)
    sample_labels = '--sample-label ' + ' '.join(sample_labels)
    library = '-l ' + library
    out = '-n mageck'
    if sgrna_len:
        sgrna_len = '--sgrna-len ' + str(sgrna_len)
    else:
        sgrna_len = ''
    if trim:
        trim = '--trim ' + str(trim)
    else:
        trim = ''

    if control_sgrna:
        control_sgrna = '--control-sgrna ' + control_sgrna
    else:
        control_sgrna = ''
    command = 'mageck count %s %s %s %s %s %s %s %s ' % (
        fastq_1, fastq_2, sample_labels, trim, library, sgrna_len, control_sgrna, out)

    logging.debug(filelist_1)
    logging.debug(filelist_2)
    logging.debug(sample_list)
    command_list = [command]

    return command_list

def mageck_count_report():
    pass
