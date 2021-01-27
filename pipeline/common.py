# -*-coding:utf-8-*-
import os
import logging


def path_check(path):
    if not os.path.exists(path):
        os.mkdir(path)
        logging.debug(path + ' created')
