__version__ = '0.0.2'

from . import preprocess as pp
from . import analysis as an
from . import plot as pl
from . import common

import logging

logging.basicConfig(filename='pipeline.log', level=logging.INFO)

# check dependency
logging.info('Checker version')

# MAGeCK
import mageck

if mageck.__version__ != '0.5.9.4':
    logging.warning('CRISPR-screen template design base on MAGeCK(0.5.9.4),your version is ' + mageck.__version__)
else:
    logging.info('MAGeCK 0.5.9.4')
