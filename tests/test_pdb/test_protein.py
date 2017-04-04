# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from dockerasmus.pdb import Protein, Chain, Atom, Residual
from dockerasmus.constants import ATOMIC_MASSES

from ..utils import DATADIR


class TestProtein(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        ## A protein made only of one arginine
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )
        ## An unrealistic protein, which properties
        ## are easy to compute manually
        cls.prot = Protein(chains={
            'A': Chain('A', residuals={
                1: Residual(1, atoms={
                    'C1': Atom(0, 0, 0, 1, 'C1'),
                    'C2': Atom(0, 0, 1, 2, 'C2')
                })
            })
        })


class TestProperties(TestProtein):

    def test_mass(self):
        self.assertEqual(self.prot.mass, 2*ATOMIC_MASSES['C'])

    def test_mass_center(self):
        self.assertEqual(list(self.prot.mass_center), [0, 0, .5])

    def test_radius(self):
        self.assertEqual(self.prot.radius, .5)


class TestMagicMethods(TestProtein):
    ## TODO: test_getitem

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

    def test_copy(self):
        arginine1 = self.arginine_prot
        arginine2 = arginine1.copy()

        for chain1, chain2 in zip(arginine1.itervalues(), arginine2.itervalues()):
            for res1, res2 in zip(chain1.itervalues(), chain2.itervalues()):
                for atom_id in res1:
                    self.assertEqual(res1[atom_id], res2[atom_id])
                    self.assertIsNot(res1[atom_id], res2[atom_id])
                    res1[atom_id].x += 1
                    self.assertNotEqual(res1[atom_id], res2[atom_id])

    def test_rmsd_ref(self):
        ## rmsd = sqrt(0.5[(0-0)²+(0-0)²+(0-0)²+(0-0)²+(0-0)²+(1-0)²])
        ##      = sqrt(0.5 * 1²) = sqrt(0.5)
        self.assertEqual(self.prot.rmsd([0,0,0]), (.5)**.5)

    def test_rmsd_prot(self):
        self.assertEqual(self.prot.rmsd(self.prot), 0)

    def test_rmsd_othertype(self):
        with self.assertRaises(TypeError):
            _ = self.prot.rmsd(1)
        with self.assertRaises(TypeError):
            _ = self.prot.rmsd("test")

    def test_contact_map_nearest(self):
        self.assertEqual( ## Nearest atom are identical
            self.prot.contact_map(self.prot, mode='nearest')[1,1],
            0,
        )

    def test_contact_map_farthest(self):
        self.assertEqual(
            self.prot.contact_map(self.prot, mode='farthest')[1,1],
            1,
        )


    def test_contact_map_mass_center(self):
        self.assertEqual( ## Same mass center
            self.prot.contact_map(self.prot, mode='mass_center')[1,1],
            0,
        )

    def test_contact_map_othermode(self):
        with self.assertRaises(ValueError):
            _ = self.prot.contact_map(self.prot, mode='rubbish')

    def test_contact_map_othertype(self):
        with self.assertRaises(TypeError):
            _ = self.prot.contact_map(1)
