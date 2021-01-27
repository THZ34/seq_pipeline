def gatk_snp(sample, kwargs):
    """"""
    reference = kwargs['reference']
    command = """gatk  MarkDuplicates -I %s.bam -O %s.markdup.bam -M %s.markdup_metrics.txt
    gatk --java-options -Xmx64G HaplotypeCaller -I %s.markup.bam -O %s.vcf -R %s""" % (
        sample, sample, sample, sample, sample, reference)
    return command


def infercnv():
    """"""
    pass
