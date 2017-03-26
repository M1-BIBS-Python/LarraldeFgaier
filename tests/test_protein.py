# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import numpy

from dockerasmus.pdb import Protein
from .utils import DATADIR


class TestProtein(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )


class TestMethods(TestProtein):

    def test_nearest_atom(self):
        pos1 = numpy.array([11.296, 86.721, 94.521]) # Coordinates of Atom 32
        self.assertEqual(self.arginine_prot.nearest_atom(pos1).id, 32)

        pos2 = numpy.array([12.924, 87.357, 96.42])
        self.assertEqual(self.arginine_prot.nearest_atom(pos2).id, 38)
