# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import numpy

from dockerasmus.pdb import Protein
from .utils import DATADIR

class TestParserOnArginine(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )
        cls.arginine = cls.arginine_prot['A'][-3]

    def test_chains(self):
        self.assertEquals(self.arginine_prot.keys(), {'A'})

    def test_residuals(self):
        self.assertEquals(self.arginine_prot['A'].keys(), {-3})

    def test_atoms(self):
        self.assertEquals(self.arginine.keys(), {
            'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD',
            'NE', 'CZ', 'NH1', 'NH2',
        })

    def test_positions(self):
        expected = {
            'C':   {'y': 86.530, 'x': 12.759, 'z': 96.365},
            'CB':  {'y': 85.746, 'x': 13.428, 'z': 93.980},
            'CA':  {'y': 85.862, 'x': 12.333, 'z': 95.041},
            'CG':  {'y': 85.172, 'x': 12.866, 'z': 92.651},
            'NE':  {'y': 85.487, 'x': 12.644, 'z': 90.195},
            'O':   {'y': 87.757, 'x': 12.924, 'z': 96.420},
            'N':   {'y': 86.721, 'x': 11.296, 'z': 94.521},
            'CZ':  {'y': 85.582, 'x': 13.114, 'z': 88.947},
            'CD':  {'y': 85.886, 'x': 13.374, 'z': 91.406},
            'NH1': {'y': 86.056, 'x': 14.338, 'z': 88.706},
            'NH2': {'y': 84.308, 'x': 14.421, 'z': 88.373},
        }

        for atom_name, pos in expected.items():
            self.assertEqual(self.arginine[atom_name].x, pos['x'])
            self.assertEqual(self.arginine[atom_name].y, pos['y'])
            self.assertEqual(self.arginine[atom_name].z, pos['z'])
            numpy.testing.assert_array_equal(
                self.arginine[atom_name].pos,
                numpy.array([pos['x'], pos['y'], pos['z']])
            )
