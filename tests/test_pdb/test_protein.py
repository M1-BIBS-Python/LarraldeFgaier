# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from dockerasmus.pdb import Protein, Chain, Atom, Residual

from ..utils import DATADIR


class TestProtein(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )


class TestMagicMethods(TestProtein):

    def test_contains_chain(self):
        self.assertIn(self.arginine_prot['A'], self.arginine_prot)
        self.assertNotIn(Chain('B'), self.arginine_prot)

    def test_contains_residual(self):
        self.assertIn(self.arginine_prot['A'][-3], self.arginine_prot)
        self.assertNotIn(Residual(-4), self.arginine_prot)

    def test_contains_atom(self):
        self.assertIn(self.arginine_prot['A'][-3]['O'], self.arginine_prot)
        self.assertNotIn(Atom(0,0,0,4), self.arginine_prot)

    def test_contains_text(self):
        self.assertIn('A', self.arginine_prot)
        self.assertNotIn('B', self.arginine_prot)

    def test_contains_othertype(self):
        with self.assertRaises(TypeError):
            _ = [1] in self.arginine_prot
        with self.assertRaises(TypeError):
            _ = set() in self.arginine_prot
        with self.assertRaises(TypeError):
            _ = object() in self.arginine_prot






class TestMethods(TestProtein):

    def test_nearest_atom(self):
        pos1 = [11.296, 86.721, 94.521] # Coordinates of Atom 32
        self.assertEqual(self.arginine_prot.nearest_atom(pos1).id, 32)

        pos2 = [12.924, 87.357, 96.42]
        self.assertEqual(self.arginine_prot.nearest_atom(pos2).id, 38)
