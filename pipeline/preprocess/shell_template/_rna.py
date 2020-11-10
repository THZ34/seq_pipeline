# Transcript quantification
def salmon(sample, kwargs):
    """Command Template:
    Transcript quantification by salmon, need quality controled .fastq file

    :param sample: sample name
    :param suffix: suffix
    :param reference: salmon reference (http://refgenomes.databio.org/)
    :return:
    """
    suffix = kwargs['suffix']
    reference = kwargs['reference']
    command = 'salmon quant ' \
              '-i %s ' \
              '-l A ' \
              '--validateMappings ' \
              '-1 fastqs/%s_qc%s ' \
              '-2 fastqs/%s_qc%s ' \
              '-o salmon/%s'% (reference, sample, suffix[0], sample, suffix[1], sample)
    return command


def htseq(samples, kwargs):
    """ Command Template:
    Transcript quantification with htseq2, need .sam file aligned by hisat2 as input

    :param samples:
    :return: command
    """
    command = 'htseq-count ' \
              '-r name ' \
              '-s no ' \
              '-a 20 ' \
              '-t exon ' \
              '-i gene_id ' \
              '-m union ' \
              '-c htseq.count ' \
              '-f sam ' \
              '-n 60 %s.sam' % ('.sam '.join(samples))
    return command

def kallisto(sample, kwargs):
    pass