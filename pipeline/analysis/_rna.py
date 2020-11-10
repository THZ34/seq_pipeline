from rpy2 import robjects
import pandas as pd
from itertools import combinations


def tximport(samples, ensttogene='ensttogene.csv', quant_type='salmon'):
    """Read salmon/kallisto result with R-tximport
    Write TPM counts to salmon/salmon_counts.csv or kallisto/kallisto_counts.csv, base on type
    Return a pandas.DataFrame contains TMP counts table
    :param samples: samples name
    :param ensttogene: The mapping table from transcript to gene symbol
    :param type: salmon or kallisto
    :return: pd.DataFrame
    """
    sample_list = []
    grp_list = []
    for sample in samples:
        if type(sample) == str:
            sample_list.append(sample)
            grp_list.append(sample)
        elif type(sample) == tuple:
            for rep in sample:
                sample_list.append(rep)
                grp_list.append(sample[0])
    r_command = """library(tximport)
    tx2gene = read.csv('%s')[c(2,3)]
    sample_list =c('%s')
    grp = c('%s')
    df = data.frame(sample_list,grp)
    fileList <- file.path("salmon", sample_list, "quant.sf")
    txi <- tximport(fileList, type = "%s", tx2gene = tx2gene)
    counts = txi$abundance
    colnames(counts) = sample_list
    write.csv(counts,'%s/%s_counts.csv')""" % (
        ensttogene, "','".join(sample_list), "','".join(grp_list), quant_type, quant_type, quant_type)
    robjects.r(r_command)
    return pd.read_csv('%s/%s_counts.csv' % (quant_type, quant_type), index_col=0)


def deg(samples, ensttogene='ensttogene.csv', quant_type='salmon'):
    sample_list = []
    grp_list = []
    for sample in samples:
        if type(sample) == str:
            sample_list.append(sample)
            grp_list.append(sample)
        elif type(sample) == tuple:
            for rep in sample:
                sample_list.append(rep)
                grp_list.append(sample[0])

    r_command_read = """library(edgeR)
                library(DESeq2)
                library(tximport)

                tx2gene = read.csv('%s')[c(2,3)]
                sample_list =c('%s')
                grp = c('%s')
                df = data.frame(sample_list,grp)
                fileList <- file.path("%s", sample_list, "quant.sf")
                txi <- tximport(fileList, type = "%s", tx2gene = tx2gene)
                dds <- DESeqDataSetFromTximport(txi, colData = df, design = ~ grp)

                d <- calcNormFactors(dds)
                cps <- cpm(d)
                k <- rowSums(cps>=1) > 2
                d <- d[k,]
                mm <- model.matrix(~-1+grp)
                d <- estimateGLMCommonDisp(d,mm)
                d <- estimateGLMTrendedDisp(d,mm)
                d <- estimateGLMTagwiseDisp(d,mm)
                f <- glmFit(d,mm)""" % (
        ensttogene, "','".join(sample_list), "','".join(grp_list), quant_type, quant_type)

    r_command_deg_list = []
    for treat, control in combinations(set(grp_list), 2):
        r_command_deg = """con <- makeContrasts("%s - %s"=grp%s-grp%s,levels=colnames(mm))
                        lrt <- glmLRT(f,contrast=con)
                        write.csv(lrt$table,'%s/%s-%s_%s.csv')""" % (
            treat, control, treat, control, treat, quant_type, control, quant_type)
        r_command_deg_list.append(r_command_deg)
    r_command = r_command_read + '\n'.join(r_command_deg_list)
    robjects.r(r_command)
    return r_command
