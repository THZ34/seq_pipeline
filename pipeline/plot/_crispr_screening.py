# -*-coding:utf-8-*-
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import math
import numpy as np


def mageck_mapped_rate(samples=None, figsize=(8, 4), save=False):
    """

    :param samples:
    :param figsize:
    :param save:
    :return:
    """
    summary = pd.read_csv('mageck.countsummary.txt', sep='\t', index_col=1)
    if not samples:
        samples = summary.index
    plt.figure(figsize=figsize)
    for sample in samples:
        total_reads = summary['Reads'][sample]
        mapped_reads = summary['Mapped'][sample]
        plt.bar(x=sample, height=total_reads, hatch='///', color='white', ec='black')
        plt.bar(x=sample, height=mapped_reads, label=sample)
        plt.text(x=sample, y=total_reads, s=str(round(mapped_reads / total_reads * 100, 1)) + '%', va='bottom',
                 ha='center', fontsize=10)
    plt.title('Mapped rates')
    plt.legend(bbox_to_anchor=(0.8, 0.5, 0.5, 0.5))
    plt.tight_layout()
    plt.xticks([])
    if save:
        plt.savefig('mageck_mapped_rate', dpi=1200)
        plt.close()
    else:
        return plt.gcf()


def mageck_abundance(samples=None, kind='Fraction', save=False, figsize=None, log_scale=False, count=None):
    """

    :param samples:
    :param kind:
    :param save:
    :param figsize:
    :param log_scale:
    :return:
    """
    if isinstance(str, count):
        count = pd.read_csv(count, sep='\t', index_col=0)
    elif isinstance(pd.DataFrame, count):
        pass
    else:
        count = pd.read_csv('mageck.count.txt', sep='\t', index_col=0)

    if not samples:
        samples = count.columns[1:]
    sample_number = len(samples)
    col_wrap = math.floor(math.sqrt(sample_number))
    if col_wrap * (col_wrap + 1) < sample_number:
        col_wrap += 1
    row_wrap = math.ceil(sample_number / col_wrap)
    if kind == 'Fraction' or kind == 'Cumulative':
        def abundance(list, kind='Cumulative'):
            list.sort()
            df = pd.DataFrame()
            cumulative = [sum(list[i:]) for i in range(len(list))]
            df['Cumulative'] = cumulative
            df['Fraction'] = df['Cumulative'] / sum(list)
            df.sort_values(by='Fraction', inplace=True, ascending=True)
            if kind == 'Cumulative':
                return df['Cumulative'].to_list()
            elif kind == 'Fraction':
                return df['Fraction'].to_list()

        if not figsize:
            figsize = (col_wrap * 4, row_wrap * 4)
        print(figsize)
        if kind == 'Fraction':
            ylabel = 'Fraction of total represented'
        else:
            ylabel = 'Cumulative count'
        if save == False:
            plt.figure(figsize=figsize)
            for i in range(sample_number):
                plt.subplot(row_wrap, col_wrap, i + 1)
                plt.plot(abundance(count.iloc[:, i + 1].to_list(), kind=kind))
                plt.title(count.columns[i + 1])
                plt.ylabel(ylabel)
                plt.xlabel('sgRNAs ranked by abundance')
                # inflection point
                x = count[count.iloc[:, i + 1] != 0].shape[0]
                if kind == 'Fraction':
                    y = 1
                elif kind == 'Cumulative':
                    y = count.iloc[i + 1].sum()
                plt.scatter(x=[x], y=[y], s=20, c='red', label='inflection point')
                plt.text(x=x, y=y * 0.95, s=str(x), va='top', ha='center')
                plt.legend(loc=2)
            plt.tight_layout()
        elif save == True:
            pdf = PdfPages('mageck_abundance.pdf')
            for sample in samples:
                plt.plot(abundance(count[sample].to_list(), kind=kind))
                plt.title(sample)
                plt.ylabel(ylabel)
                plt.xlabel('sgRNAs ranked by abundance')
                # inflection point
                x = count[count[sample] != 0].shape[0]
                if kind == 'Fraction':
                    y = 1
                elif kind == 'Cumulative':
                    y = count[sample].sum()
                plt.scatter(x=[x], y=[y], s=20, c='red', label='inflection point')
                plt.text(x=x, y=y * 0.95, s=str(x), va='top', ha='center')
                plt.legend(loc=2)
                pdf.savefig()
                plt.close()
            pdf.close()
            plt.figure(figsize=figsize)
            for i in range(sample_number):
                plt.subplot(row_wrap, col_wrap, i + 1)
                plt.plot(abundance(count.iloc[:, i + 1].to_list(), kind=kind))
                plt.title(count.columns[i + 1])
                plt.ylabel(ylabel)
                plt.xlabel('sgRNAs ranked by abundance')
                # inflection point
                x = count[count.iloc[:, i + 1] != 0].shape[0]
                if kind == 'Fraction':
                    y = 1
                elif kind == 'Cumulative':
                    y = count.iloc[i + 1].sum()
                plt.scatter(x=[x], y=[y], s=20, c='red', label='inflection point')
                plt.text(x=x, y=y * 0.95, s=str(x), va='top', ha='center')
                plt.legend(loc=2)
            plt.tight_layout()
            plt.savefig('mageck_abundance', dpi=1200)
            plt.close()
    elif kind == 'bar':
        if not figsize:
            figsize = (8, 4)
        plt.figure(figsize=figsize)
        for sample in samples:
            plt.bar(x=sample, height=count[count[sample] != 0].shape[0], label=sample)
            plt.text(x=sample, y=count[count[sample] != 0].shape[0], s=str(count[count[sample] != 0].shape[0]),
                     ha='center', va='bottom')
        plt.xticks([])
        plt.legend(bbox_to_anchor=(0.8, 0.5, 0.5, 0.5))
        plt.title('Number of sgRNA detected')
        plt.tight_layout()
        if save == True:
            plt.savefig('mageck_detected_sgrna', dpi=1200)
        plt.close()
    elif kind == 'dist':
        restruct_count = pd.DataFrame()
        for sample in count.columns[1:]:
            tmp_df = count[['Gene', sample]]
            tmp_df['sample'] = sample
            tmp_df.columns = ['Gene', 'sgRNA count', 'sample']
            restruct_count = pd.concat([restruct_count, tmp_df], axis=0, sort=False)
        if log_scale:
            restruct_count['sgRNA count'] = restruct_count['sgRNA count'] + 1
        if save == False:
            g = sns.FacetGrid(restruct_count, col='sample', col_wrap=col_wrap)
            g.map(sns.histplot, 'sgRNA count', log_scale=log_scale, bins=40)
        elif save == True:
            g = sns.FacetGrid(restruct_count, col='sample', col_wrap=col_wrap)
            g.map(sns.histplot, 'sgRNA count', log_scale=log_scale, bins=40)
            plt.savefig('mageck_distribution', dpi=1200)
            plt.close()
            pdf = PdfPages('mageck_distribution.pdf')
            for sample in samples:
                sns.histplot(restruct_count[restruct_count['sample'] == sample], log_scale=log_scale, bins=40)
                pdf.savefig()
                plt.close()
            pdf.close()


def mageck_saturation(samples=None, save=False, sampling=None, seed=1, count=None):
    """

    :return:
    """
    if isinstance(str, count):
        count = pd.read_csv(count, sep='\t', index_col=0)
    elif isinstance(pd.DataFrame, count):
        pass
    else:
        count = pd.read_csv('mageck.count.txt', sep='\t', index_col=0)
    count_weight = count.iloc[:, 1:] / count.iloc[:, 1:].sum()
    if not samples:
        samples = count.columns[1:]
    sample_number = len(samples)
    col_wrap = math.floor(math.sqrt(sample_number))
    if col_wrap * (col_wrap + 1) < sample_number:
        col_wrap += 1
    if not sampling:
        sampling = range(0, 200000, 10000)
    saturation_df = pd.DataFrame()
    np.random.seed(seed)
    for sample in samples:
        counts_sgrnas = []
        for choice_number in sampling:
            sub_samples = np.random.choice(count_weight.index, choice_number, p=count_weight[sample])
            counts_sgrnas.append(len(set(sub_samples)))
        tmp_df = pd.DataFrame([counts_sgrnas, [sample] * len(sampling), list(sampling)],
                              index=['Detected sgRNA number', 'sample', 'Count number']).T
        saturation_df = pd.concat([saturation_df, tmp_df], axis=0, sort=False)
    saturation_df['Detected sgRNA number'] = saturation_df['Detected sgRNA number'].astype(dtype=int)
    saturation_df['Count number'] = saturation_df['Count number'].astype(dtype=int)
    if save == True:
        pdf = PdfPages('mageck_saturation.pdf')
        sns.lineplot(data=saturation_df, x='Count number', y='Detected sgRNA number', hue='sample', marker='^')
        pdf.savefig()
        plt.savefig('mageck_saturation', dpi=1200)
        plt.close()
        g = sns.FacetGrid(data=saturation_df, col='sample', col_wrap=col_wrap)
        g.map(sns.lineplot, 'Count number', 'Detected sgRNA number')
        pdf.savefig()
        plt.close()
        for sample in samples:
            sns.lineplot(data=saturation_df[saturation_df['sample'] == sample], x='Count number',
                         y='Detected sgRNA number', hue='sample', marker='^')
            pdf.savefig()
            plt.close()
        pdf.close()
        saturation_df.to_csv('mageck_saturation.csv')
    elif save == False:
        sns.lineplot(data=saturation_df, x='Count number', y='Detected sgRNA number', hue='sample', marker='^')
