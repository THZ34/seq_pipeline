from rpy2 import robjects
import pandas as pd


def tximport(samples, ):
    """"""

    sample_list = []
    grp_list = []
    for sample in samples:
        if type(sample)==str:
            sample_list.append(sample)
            grp_list.append(sample)
        elif type(sample)==tuple:
            for rep in sample:
                sample_list.append(rep)
                grp_list.append(sample[0])
    r_command = """tx2gene = read.csv('ensttogene.csv')[c(2,3)]
                    sample_list =c('%s')
                    grp = c('%s')
                    df = data.frame(sample_list,grp)
                    fileList <- file.path("salmon", sample_list, "quant.sf")
                    txi <- tximport(fileList, type = "salmon", tx2gene = tx2gene)
                    counts = txi$counts
                    colnames(counts) = sample_list
                    write.csv(counts,'salmon_counts.csv')"""%()

