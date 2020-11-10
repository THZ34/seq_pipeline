
# macs2 cell peak
def macs2_cuttag(sample,kwargs):
    """Shell command template
    Call peak with MACS2
    This template for cut&tag

    :param sample:
    :return: command
    """
    command = 'macs2 callpeak ' \
              '-t sam/%s.sam ' \
              '--keep-dup all ' \
              '-n macs2/%s ' \
              '-B ' \
              '-p 1e-5 ' \
              '-g hs' % (sample, sample)

    return command

