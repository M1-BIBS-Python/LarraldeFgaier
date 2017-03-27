# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from dockerasmus.score.forcefield import cornell
from dockerasmus.pdb import Protein

from ..utils import DATADIR


class TestCornellScoringFunction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.barnase = Protein.from_pdb_file(os.path.join(DATADIR, "barnase.native.pdb.gz"))
        cls.barstar = Protein.from_pdb_file(os.path.join(DATADIR, "barstar.native.pdb.gz"))

    @unittest.skipIf(cornell.score._backend=="numpy", "theano is not available")
    def test_backend_agnostic(self):
        score_t = cornell.score                            # precompiled theano
        score_n = cornell._CornellScoringFunction('numpy') # numpy backend

        self.assertAlmostEqual(
            score_t(self.barnase, self.barstar),
            score_n(self.barnase, self.barstar),
            places=2,
        )
