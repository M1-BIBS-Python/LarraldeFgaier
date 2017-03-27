# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from dockerasmus.pdb import Protein, Residual, Atom

from .utils import DATADIR


class TestResidual(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )['A'][-3]


class TestProperties(TestResidual):

    def test_mass_center(self):
        # Actual mass center were calculated manually using AARG
        # as the reference position of the atoms.
        self.assertAlmostEqual(self.arginine.mass_center[0], 12.7554067, 6)
        self.assertAlmostEqual(self.arginine.mass_center[1], 85.3782889, 6)
        self.assertAlmostEqual(self.arginine.mass_center[2], 91.9005718, 6)


class TestContains(TestResidual):

    def test_contains_id(self):
        self.assertIn(32, self.arginine)
        self.assertNotIn(75, self.arginine)

    def test_contains_element(self):
        nitrogen = self.arginine['N']
        oxygen = Atom(13.753, 87.471, 95.27, 38, 'O')
        sulfur = Atom(0, 0, 0, 1, "S")
        self.assertIn(nitrogen, self.arginine)
        self.assertIn(oxygen, self.arginine)
        self.assertNotIn(sulfur, self.arginine)

    def test_contains_name(self):
        self.assertIn('N', self.arginine)
        self.assertIn('CZ', self.arginine)
        self.assertNotIn('SNCF', self.arginine)

    def test_contains_error(self):
        with self.assertRaises(TypeError):
            [] in self.arginine
