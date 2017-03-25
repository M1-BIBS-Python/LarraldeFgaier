# coding: utf-8
"""
Test doctest contained tests in every file of the module.
"""
from __future__ import (
    absolute_import,
    unicode_literals,
)

import os
import sys
import doctest
import re
import warnings

import dockerasmus


class IgnoreUnicodeChecker(doctest.OutputChecker):
    def check_output(self, want, got, optionflags):
        if sys.version_info[0] > 2:
            want = re.sub("u'(.*?)'", "'\\1'", want)
            want = re.sub('u"(.*?)"', '"\\1"', want)
        return doctest.OutputChecker.check_output(self, want, got, optionflags)


def _load_tests_from_module(tests, module, globs, setUp, tearDown, recurse=5):
    """Load tests from module, iterating through submodules"""
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith('_')):
        if isinstance(attr, type(os)):
            if attr.__name__.startswith(module.__name__.split('.')[0]):
                try:
                    tests.addTests(doctest.DocTestSuite(attr, globs=globs,
                        setUp=setUp, tearDown=tearDown, checker=IgnoreUnicodeChecker()))
                except ValueError:
                    pass
                else:
                    if recurse:
                        tests=_load_tests_from_module(tests, attr, globs, setUp, tearDown, recurse-1)
    return tests


def load_tests(loader, tests, ignore):
    """load_test function used by unittest to find the doctests"""

    def _setUp(self):
        """setUp method used by the DocTestSuite"""
        pass

    def _tearDown(self):
        """tearDown method used by the DocTestSuite"""
        pass

    globs = {
        # globs for dockerasmus
        'dockerasmus': dockerasmus,

        # globs for dockerasmus.utils
        'nth': dockerasmus.utils.iterators.nth,
        'infinitewords': dockerasmus.utils.iterators.infinitewords,
        'method_requires': dockerasmus.utils.decorators.method_requires,

        # globs for pdb:
        'Protein': dockerasmus.pdb.Protein,

    }
    tests = _load_tests_from_module(tests, dockerasmus, globs, _setUp, _tearDown)
    return tests


def setUpModule():
    warnings.simplefilter('ignore')


def tearDownModule():
    warnings.simplefilter(warnings.defaultaction)
