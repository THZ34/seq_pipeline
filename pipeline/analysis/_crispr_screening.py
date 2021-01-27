from pipeline.common import path_check
import pandas as pd
import os
import logging
import multiprocessing
import math
import numpy as np

cpu = math.floor(multiprocessing.cpu_count() * 0.8)


def mageck_rra(sample_pairs=None):
    """Shell command template
    Run Mageck test

    :param sample_pairs: (control,treatment)
    :return:
    """
    # check directory
    path_check('mageck_rramkdi')

    commands = []
    if not sample_pairs:
        command = 'mageck test ' \
                  '-k mageck.count_normalized.txt ' \
                  '--paired ' \
                  '--pdf-report ' \
                  '-n mageck/paired '

        commands.append(command)
    else:
        for pair in sample_pairs:
            command = 'mageck test ' \
                      '-k mageck.count_normalized.txt ' \
                      '-t %s ' \
                      '-c %s ' \
                      '--pdf-report ' \
                      '--norm-method none ' \
                      '-n mageck_rra/%s-vs-%s ' % (pair[1], pair[0], pair[1], pair[0])
            commands.append(command)
    return commands


def mageck_mle(design_matrix):
    """

    :param design_matrix:
    :return:
    """
    # check directory
    path_check('mageck_mle')

    # args parse
    if isinstance(design_matrix, pd.DataFrame):
        pass
    elif isinstance(design_matrix, str):
        os.path.exists(design_matrix)
        design_matrix = pd.read('design_matrix.xlsx', index_col=0)
    else:
        logging.error('Design matrix must be pd.DataFrame or EXCEL file')
    commands = []
    for control in design_matrix.index:
        path_check('mageck_mle/' + control)
        samples = design_matrix.loc[control][design_matrix.loc[control] != 0].index
        mageck_design = pd.DataFrame(np.eye(len(samples) + 1))
        columns = ['baseline']
        columns.extend(samples)
        index = [control]
        index.extend(samples)
        mageck_design.columns = columns
        mageck_design.index = index
        mageck_design.index.name = 'Samples'
        mageck_design['baseline'] = 1
        mageck_design = mageck_design.astype(dtype=int)
        mageck_design.to_csv('mageck_mle/%s/mageck_mle_design.txt' % control, sep='\t')
        command = 'mageck mle -k mageck.count_normalized.txt -d mageck_mle/%s/mageck_mle_design.txt -n mageck_mle/%s/%s --norm-method none --threads %s' % (
            control, control, control, cpu)
        commands.append(command)
    return commands
