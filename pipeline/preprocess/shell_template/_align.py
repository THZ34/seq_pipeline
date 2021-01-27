import multiprocessing
import math


def bwa(samples, kwargs):
    """

    :param sample:
    :param kwargs:
    :return:
    """
    suffix = kwargs['suffix']
    reference = kwargs['reference']
    command_list = []
    for sample in samples:
        if type(sample) == str:
            command = """bwa mem -t %s -R '@RG\tID:%s\tPL:illumina\tSM:%s' %s %s_qc%s %s_qc%s | samtools view -b |samtools sort > sam/%s.bam
            samtools index sam/%s.bam""" % (
                max(1, math.floor(multiprocessing.cpu_count() * 0.9)), sample, sample, reference, sample, suffix[0],
                sample, suffix[1], sample, sample)
            command_list.append(command)
        elif type(sample) == tuple:
            for rep in sample:
                command = """bwa mem -t %s -R '@RG\tID:%s\tPL:illumina\tSM:%s' %s %s_qc%s %s_qc%s | samtools view -b | samtools sort > sam/%s.bam
                samtools index sam/%s.bam""" % (
                    max(1, math.floor(multiprocessing.cpu_count() * 0.9)), rep, rep, reference, sample[0], suffix[0],
                    rep, suffix[1], rep, rep)
                command_list.append(command)

    return command_list


def bowtie2_cuttag(sample, kwargs):
    """Alignment with bowtie2
    This template for cut&tag

    :param sample: sample name
    :param suffix: suffix
    :param reference: bowtie2 reference (http://refgenomes.databio.org/)
    :return: command
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
              '-1 fastqs/%s_qc%s ' \
              '-2 fastqs/%s_qc%s ' \
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
    :param sample: sample name
    :return: command
    """
    command = 'samtools view -b sam/%s.sam | samtools sort > sam/%s.sorted.bam \n' \
              'samtools index sam/%s.sorted.bam' % (sample, sample, sample)
    return command


# bam coverage : bigwig
def bam_coverage(sample, kwargs):
    """Shell command template
    Calculate reads coverage

    :param sample:
    :return: command
    """
    command = 'bamCoverage -b sam/%s.sorted.bam -o sam/%s.bigwig -of bigwig' % (sample, sample)
    return command


# hisat2 RNA-seq
def hisat2(sample, kwargs):
    """Shell command template
    Alignment with hisat2

    :param sample: sample name
    :param suffix: suffix
    :param reference: hisat2 reference
    :return: command
    """
    suffix = kwargs['suffix']
    reference = kwargs['reference']
    command = 'hisat2 ' \
              '-p 60 ' \
              '-x %s ' \
              '-1 fastqs/%s_qc%s ' \
              '-2 fastqs/%s_qc%s ' \
              '-S sam/%s.sam' % (reference, sample, suffix[0], sample, suffix[1], sample)
    return command

# STAR alignment
def star(sample, kwargs):
    """Shell command template
     Alignment with star

     :param samples: sample name
     :param suffix: suffix
     :param reference: star reference
     :return: command
     """
    suffix = kwargs['suffix']
    reference=kwargs['reference']
    command = 'STAR ' \
              '--runThreadN 60 ' \
              '--genomeDir %s ' \
              '--readFilesIn fastqs/%s_qc%s fastqs/%s_qc%s '\
              '--outSAMtype BAM SortedByCoordinate '\
              '--outFileNamePrefix star/%s' % (reference, sample, suffix[0], sample, suffix[1], sample)
    return command
