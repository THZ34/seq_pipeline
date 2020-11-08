
# macs2 cell peak
def macs2(sample,kwargs):
    """Shell command template
    Call peak with MACS2

    :param sample:
    :return:
    """
    command = 'macs2 callpeak ' \
              '-t sam/%s.sam ' \
              '--keep-dup all ' \
              '-n macs2/%s ' \
              '-B ' \
              '-p 1e-5 ' \
              '-g hs' % (sample, sample)

    return command

