# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest
import numpy

from dockerasmus.pdb import Protein

from ..utils import DATADIR


class TestParserOnArginine(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )
        cls.arginine = cls.arginine_prot['A'][-3]

    def test_chains(self):
        self.assertEqual(set(self.arginine_prot.keys()), {'A'})

    def test_residuals(self):
        self.assertEqual(set(self.arginine_prot['A'].keys()), {-3})

    def test_atoms(self):
        self.assertEqual(
            set(self.arginine.keys()),
            {'N', 'CA', 'C', 'O', 'CB', 'CG',
             'CD', 'NE', 'CZ', 'NH1', 'NH2'}
        )

    def test_positions(self):
        expected = {
            'N':   {'x': 11.281, 'y': 86.699, 'z': 94.383},
            'CA':  {'x': 12.353, 'y': 85.696, 'z': 94.456},
            'C':   {'x': 13.559, 'y': 86.257, 'z': 95.222},
            'O':   {'x': 13.753, 'y': 87.471, 'z': 95.270},
            'CB':  {'x': 12.774, 'y': 85.306, 'z': 93.039},
            'CG':  {'x': 11.754, 'y': 84.432, 'z': 92.321},
            'CD':  {'x': 11.698, 'y': 84.678, 'z': 90.815},
            'NE':  {'x': 12.984, 'y': 84.447, 'z': 90.163},
            'CZ':  {'x': 13.202, 'y': 84.534, 'z': 88.850},
            'NH1': {'x': 12.218, 'y': 84.840, 'z': 88.007},
            'NH2': {'x': 14.421, 'y': 84.308, 'z': 88.373},
        }



        for atom_name, pos in expected.items():
            self.assertEqual(self.arginine[atom_name].x, pos['x'])
            self.assertEqual(self.arginine[atom_name].y, pos['y'])
            self.assertEqual(self.arginine[atom_name].z, pos['z'])
            numpy.testing.assert_array_equal(
                self.arginine[atom_name].pos,
                numpy.array([pos['x'], pos['y'], pos['z']])
            )
