def bwa(samples, suffix, reference):
    pass


def bowtie2_cuttag(sample, kwargs):
    """

    :param sample:
    :param suffix:
    :param reference:
    :return:
    """
    suffix = kwargs['suffix']
    reference = kwargs['reference']
    command = 'bowtie2 ' \
              '-p 60 ' \
              '--local ' \
              '--very-sensitive-local ' \
              '--no-unal ' \
              '--no-mixed ' \
              '--no-discordant ' \
              '--phred33 ' \
              '-I 10 ' \
              '-X 700 ' \
              '-x %s ' \
              '-1 fastq/%s_qc%s ' \
              '-2 fastq/%s_qc%s ' \
              '-S sam/%s.sam' % (reference, sample, suffix[0], sample, suffix[1], sample)
    return command


def bowtie2_chip(samples, suffix, reference):
    """"""
    pass


def bowtie2_atac(samples, suffix, reference):
    pass


# sort and index
def sam_sort(sample, kwargs):
    """Shell command template
    Convert sam to bam, then sort and make index
    :param sample:
    :param kwargs:
    :return:
    """
    command = 'samtools view -b sam/%s.sam | samtools sort > sam/%s.sorted.bam \n' \
              'samtools index sam/%s.sorted.bam' % (sample, sample, sample)
    return command


# bam coverage : bigwig
def bam_coverage(sample, kwargs):
    """Shell command template
    Calculate reads coverage

    :param sample:
    :param kwargs:
    :return:
    """
    command = 'bamCoverage -b sam/%s.sorted.bam -o sam/%s.bigwig -of bigwig' % (sample, sample)
    return command


# hisat2 RNA-seq
def hisat2(sample, kwargs):
    """Shell command template
    Alignment with hisat2

    :param sample:
    :param reference:
    :return:
    """
    suffix = kwargs['suffix']
    reference = kwargs['reference']
    command = 'hisat2 ' \
              '-p 60 ' \
              '-x %s ' \
              '-1 fastq/%s_qc%s ' \
              '-2 fastq/%s_qc%s ' \
              '-S sam/%s.sam' % (reference, sample, suffix[0], sample, suffix[1], sample)
    return command
