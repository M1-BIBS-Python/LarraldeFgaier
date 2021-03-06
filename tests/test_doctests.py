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
import numpy
import theano
import types

import dockerasmus
import dockerasmus.score
import dockerasmus.pdb


from .utils import DATADIR


class IgnoreUnicodeChecker(doctest.OutputChecker):
    """A checker that removes the 'u' string prefix from expected results
    """
    def check_output(self, want, got, optionflags):
        if sys.version_info[0] > 2:
            want = re.sub("u'(.*?)'", "'\\1'", want)
            want = re.sub('u"(.*?)"', '"\\1"', want)
        return doctest.OutputChecker.check_output(self, want, got, optionflags)


def _load_tests_from_module(tests, module, globs, setUp, tearDown, recurse=5):
    """Load tests from module, iterating through submodules"""
    for attr in (getattr(module, x) for x in dir(module) if not x.startswith('_')):
        if isinstance(attr, types.ModuleType):
            if attr.__name__.startswith(module.__name__.split('.')[0]):
                try:
                    tests.addTests(doctest.DocTestSuite(
                        attr, globs=globs,
                        optionflags=doctest.ELLIPSIS, setUp=setUp,
                        tearDown=tearDown, checker=IgnoreUnicodeChecker()
                    ))
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
        numpy.set_printoptions(precision=3)

    def _tearDown(self):
        """tearDown method used by the DocTestSuite"""
        pass

    globs = {
        # generic modules
        'dockerasmus': dockerasmus,
        'numpy': numpy,
        'theano': theano,

        # globs for dockerasmus.utils
        'nth': dockerasmus.utils.iterators.nth,
        'wordrange': dockerasmus.utils.iterators.wordrange,
        'method_requires': dockerasmus.utils.decorators.method_requires,
        'distance': dockerasmus.utils.matrices.distance,
        'normalized': dockerasmus.utils.matrices.normalized,
        'maybe_import': dockerasmus.utils.maybe_import,

        # globs for dockerasmus.score
        'ScoringFunction': dockerasmus.score.ScoringFunction,
        'LennardJones': dockerasmus.score.components.LennardJones,
        'Fabiola': dockerasmus.score.components.Fabiola,
        'Coulomb': dockerasmus.score.components.Coulomb,

        # globs for pdb:
        'Protein': dockerasmus.pdb.Protein,

        # locals
        'arginine': dockerasmus.pdb.Protein.from_pdb_file(os.path.join(DATADIR, 'arginine.pdb'))['A'][-3],
        'barnase': dockerasmus.pdb.Protein.from_pdb_file(os.path.join(DATADIR, 'barnase.native.pdb.gz')),
        'barstar': dockerasmus.pdb.Protein.from_pdb_file(os.path.join(DATADIR, 'barstar.native.pdb.gz'))
    }
    tests = _load_tests_from_module(tests, dockerasmus, globs, _setUp, _tearDown)
    return tests


def setUpModule():
    warnings.simplefilter('ignore')


def tearDownModule():
    warnings.simplefilter(warnings.defaultaction)
