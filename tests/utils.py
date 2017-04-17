# coding: utf-8
from __future__ import unicode_literals

import os
import functools

CURRDIR = os.path.dirname(os.path.abspath(__file__))
DATADIR = os.path.join(CURRDIR, "data")

try:
    from unittest import mock
except ImportError: # pragma: no cover
    import mock


def suppress_tf_log(func):
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        with mock.patch.dict('os.environ', {'TF_CPP_MIN_LOG_LEVEL': "3"}):
            return func(*args, **kwargs)
    return new_func
