# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

__author__ = "Martin Larralde (althonos)"
__author_email__ = "martin.larralde@ens-cachan.fr"
__version__ = "0.1.1"
__license__ = "GPLv3"

try:
    from . import pdb
    from . import score
    from . import utils
except ImportError as import_error:
    import warnings
    warnings.warn(str(import_error), ImportWarning)
